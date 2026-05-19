#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <memory>
#include <map>
#include <sstream>
#include <array>

// Assuming these types and constants are defined in other headers
// For the purpose of this translation, we provide stubs or assume they exist.
// In a real project, these would be in separate header files.

// Constants
constexpr int BUFSIZE = 256;
constexpr int RKIND = 8; // Assuming double precision
using REAL = double;
using INT4 = int32_t;
using INT8 = int64_t;

// Mock types to make the code compile structurally
// These correspond to Fortran derived types used in the module

struct XYZlimit_t {
    INT4 XI, XE, YI, YE, ZI, ZE;
};

struct xyzlimit_scaled_t {
    XYZlimit_t base;
    REAL scale;
};

struct limit_t {
    XYZlimit_t bounds;
    // Other fields omitted for brevity as they are not directly accessed in the snippet
};

struct tagtype_t {
    // Placeholder
};

struct taglist_t {
    struct {
        std::vector<std::vector<std::vector<INT4>>> x;
        std::vector<std::vector<std::vector<INT4>>> y;
        std::vector<std::vector<std::vector<INT4>>> z;
    } edge;
    struct {
        std::vector<std::vector<std::vector<INT4>>> x;
        std::vector<std::vector<std::vector<INT4>>> y;
        std::vector<std::vector<std::vector<INT4>>> z;
    } face;
};

struct MedioExtra_t {
    bool exists = false;
    INT4 index = 0;
};

struct PlaneWave_t {
    bool isRC = false;
    INT4 numModes = 0;
    REAL incertmax = 0.0;
    std::vector<REAL> px, py, pz, ex, ey, ez, INCERT;
};

struct SGGFDTDINFO_t {
    bool thereAreMagneticMedia = false;
    bool thereArePMLMagneticMedia = false;
    REAL dt = 0.0;
    REAL DX = 0.0, DY = 0.0, DZ = 0.0;
    
    struct {
        INT4 Conta = 0;
        INT4 MaxConta = 0;
        std::vector<int> elem;
    } EShared;

    INT4 NumMedia = 0;
    INT4 AllocMed = 0;
    // Assuming Med is a vector of media objects
    std::vector<int> Med; 

    struct {
        INT4 XI, XE, YI, YE, ZI, ZE;
    } Alloc[6]; // Simplified access to Alloc array

    INT4 NumPlaneWaves = 0;
    std::vector<PlaneWave_t> PlaneWave;
};

struct media_matrices_t {
    // 3D arrays for tags and indices
    // Using flattened vectors or vectors of vectors for 3D
    std::vector<std::vector<std::vector<INT4>>> sggMtag;
    std::vector<std::vector<std::vector<INT4>>> sggMiNo;
    std::vector<std::vector<std::vector<INT4>>> sggMiEx;
    std::vector<std::vector<std::vector<INT4>>> sggMiEy;
    std::vector<std::vector<std::vector<INT4>>> sggMiEz;
    std::vector<std::vector<std::vector<INT4>>> sggMiHx;
    std::vector<std::vector<std::vector<INT4>>> sggMiHy;
    std::vector<std::vector<std::vector<INT4>>> sggMiHz;
};

// Mock Parser/Preprocessor class components
struct ConformalRegs_t {
    bool hasVolumes = false;
};

struct PlnSrc_t {
    INT4 nc = 0;
    struct {
        REAL coor1[3];
        REAL coor2[3];
        bool isRC = false;
        INT4 numModes = 0;
        REAL incertmax = 0.0;
    } collection[10]; // Placeholder size
};

struct Parseador_t {
    ConformalRegs_t conformalRegs;
    PlnSrc_t plnSrc;
    // Other regions omitted for brevity
    struct { INT4 nvols=0, nsurfs=0, nLINS=0; } pmcregs, DielRegs, FRQDEPMATS, ANIMATS, LossyThinSurfs, twires, swires, tSlots;
};

// Helper functions stubs
void print11(INT4 layoutnumber, const std::string& msg) {
    std::cout << "[Proc " << layoutnumber << "] " << msg << std::endl;
}

void cuentatags(Parseador_t& this_, tagtype_t& tagtype, INT4 layoutnumber, const std::string& fichin) {
    /* TODO: Count tags from input file. From preprocess_geom.F90 */
}

void prepro_skindepth(Parseador_t& this_, const std::string& fichin) {
    /* TODO: Preprocess skin depth from input file. From preprocess_geom.F90 */
}

void populatePlaneWaveRC(PlaneWave_t& pw, bool simu_devia) {
    /* TODO: Populate plane wave random constants. From preprocess_geom.F90 */
}

// MPI stubs
#ifdef CompileWithMPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 0
#define MPI_SUM 0
#define MPI_BARRIER 0
#define MPI_AllReduce(a,b,c,d,e,f,g) do{}while(0)
#define MPI_Barrier(a,b) do{}while(0)
#endif

// Conformal helpers stubs
struct ConformalMedia_t {
    REAL time_step_scale_factor = 1.0;
};

std::vector<ConformalMedia_t> buildConformalMedia(const ConformalRegs_t& regs) {
    return {};
}

struct side_tris_map_t {};
std::vector<side_tris_map_t> buildSideMaps(const ConformalRegs_t& regs) {
    return {};
}

std::vector<REAL> getDifferentEdgeRatios(const std::vector<ConformalMedia_t>& cm) {
    return {1.0};
}

std::vector<REAL> getDifferentFaceRatios(const std::vector<ConformalMedia_t>& cm) {
    return {1.0};
}

// Constants for field indices (mock)
enum FieldIndex {
    iHx, iHy, iHz, iEx, iEy, iEz
};

namespace Preprocess_m {

    // Global variables
    REAL cluz = 0.0;
    REAL zvac = 0.0;
    REAL eps0 = 0.0;
    REAL mu0 = 0.0;

    // Function declarations
    void read_geomData(
        media_matrices_t& media,
        const std::string& fichin,
        INT4 layoutnumber,
        INT4 num_procs,
        std::vector<limit_t>& SINPML_fullsize,
        std::vector<limit_t>& fullsize,
        SGGFDTDINFO_t& sgg,
        std::vector<limit_t>& groundwires,
        REAL attfactor,
        bool& mibc,
        bool& SGBC,
        bool& SGBCDispersive,
        MedioExtra_t& MEDIOEXTRA,
        REAL maxSourceValue,
        REAL skindepthpre,
        bool& createmapvtk,
        bool input_conformal_flag,
        bool& CLIPREGION,
        REAL boundwireradius,
        REAL maxwireradius,
        bool updateshared,
        bool run_with_dmma,
        REAL eps00,
        REAL mu00,
        bool simu_devia,
        bool hay_slanted_wires,
        bool verbose,
        bool ignoresamplingerrors,
        tagtype_t& tagtype,
        const std::string& wiresflavor
    );

    void read_limits_nogeom();
    void AssigLossyOrPECtoNodes();
    void searchtag();
    void checkDielectricTagForDuplicate();
    void checkAnimatedTagForDuplicate();
    void checkLossyTagForDuplicate();

    // Implementation
    void read_geomData(
        media_matrices_t& media,
        const std::string& fichin,
        INT4 layoutnumber,
        INT4 num_procs,
        std::vector<limit_t>& SINPML_fullsize,
        std::vector<limit_t>& fullsize,
        SGGFDTDINFO_t& sgg,
        std::vector<limit_t>& groundwires,
        REAL attfactor,
        bool& mibc,
        bool& SGBC,
        bool& SGBCDispersive,
        MedioExtra_t& MEDIOEXTRA,
        REAL maxSourceValue,
        REAL skindepthpre,
        bool& createmapvtk,
        bool input_conformal_flag,
        bool& CLIPREGION,
        REAL boundwireradius,
        REAL maxwireradius,
        bool updateshared,
        bool run_with_dmma,
        REAL eps00,
        REAL mu00,
        bool simu_devia,
        bool hay_slanted_wires,
        bool verbose,
        bool ignoresamplingerrors,
        tagtype_t& tagtype,
        const std::string& wiresflavor
    ) {
        // Local variables
        INT4 tama, tama2, tama3, tama4, tama5, tama6, i, j, k, tipotemp, tamaSonda, tamaoldSONDA, tamaBloquePrb, tamaScrPrb, pozi, tama2bis, numeroasignaciones, ci;
        INT4 orientacion, orientacionL, orientacionR, direccion, contamedia, oldcontamedia, maxcontamedia, mincontamedia, inicontamedia, i1, j1, field, k1, pecmedio, ii, medio1, medio2, sondas, CONTACURR, CONTAVOLT, I_, J_;
        REAL ex, ey, ez, px, py, pz, amplitud, minSpaceStep;
        XYZlimit_t punto, BoundingBox, conf_bounding_box;
        INT4 nsurfs, numus, OrigIndex, numminus;
        REAL delta, del, sig_max;
        INT4 conta1, conta2, MEDIO, imenos1, jmenos1, kmenos1, o, p, puntoxi, puntoyi, puntozi;
        INT4 bboxwirXI, dummy_bboxwirXI, bboxwirYI, dummy_bboxwirYI, bboxwirzI, dummy_bboxwirzI;
        INT4 bboxwirXE, dummy_bboxwirXE, bboxwirYE, dummy_bboxwirYE, bboxwirZE, dummy_bboxwirZE, IERR;
        INT8 memo;
        
        std::string whoami, whoamishort, ext, extpoint, chari, charj, chark, chari2, charj2, chark2, MultiportFile, MultiportFile2, buff;
        std::vector<REAL> dummy_px, dummy_py, dummy_pz, dummy_ex, dummy_ey, dummy_ez, dummy_INCERT;
        
        bool isathinwire, VALIDO, existia, medioespecial, nodo_cazado, errnofile, errnofile1, errnofile2, errnofile3, errnofile4;
        REAL tiempo1, tiempo2, field1, field2, rdummy;
        std::vector<INT4> contapuntos;
        
        bool paraerrhilo, groundwires_flag, islossy, DENTRO;
        REAL width, dir[3], epr1, mur1;
        bool oriX, oriY, oriZ, oriX2, oriY2, oriZ2, oriX3, oriY3, oriZ3, iguales;
        bool oriX4, oriY4, oriZ4;
        REAL EprSlot[3][3], MurSlot[3][3];
        INT4 indicemedio, i11, j11;
        
        INT4 numertag;
        INT4 Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE;
        INT4 Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE;
        INT4 Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE;
        INT4 Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE;
        INT4 Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE;
        INT4 Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE;

        // Initialize globals
        eps0 = eps00;
        mu0 = mu00;
        cluz = 1.0 / std::sqrt(eps0 * mu0);
        zvac = std::sqrt(mu0 / eps0);

        // Call cuentatags
        // Note: In C++, we might need to pass references or pointers depending on implementation
        // Assuming tagtype is passed by reference or modified internally
        tagtype_t local_tagtype;
        cuentatags(/*this*/ nullptr, local_tagtype, layoutnumber, fichin); // Placeholder for 'this'

        delta = -1.0;

        sgg.thereAreMagneticMedia = true;
        sgg.thereArePMLMagneticMedia = true;

        if (skindepthpre) {
            if (layoutnumber == 0) {
                print11(layoutnumber, "Preprocessing SGBC materials to include skin-depth effects....");
                prepro_skindepth(/*this*/ nullptr, fichin); // Placeholder
                print11(layoutnumber, "Finished preprocessing for skin-depth.");
            }
#ifdef CompileWithMPI
            MPI_Barrier(MPI_COMM_WORLD, &IERR);
#endif
            return;
        }

        // Create whoami string
        char whoami_buf[BUFSIZE];
        snprintf(whoami_buf, BUFSIZE, "(%d/%d) ", layoutnumber + 1, num_procs);
        whoami = whoami_buf;
        
        char whoamishort_buf[BUFSIZE];
        snprintf(whoamishort_buf, BUFSIZE, "%d", layoutnumber + 1);
        whoamishort = whoamishort_buf;

        // Create space for etangential shared info
        sgg.EShared.Conta = 0;
        sgg.EShared.MaxConta = 10;
        sgg.EShared.elem.resize(10);

        // Conformal Media Block
        {
            INT4 m;
            REAL min_scale_factor = 1.0, dt;
            // Assuming 'this' has conformalRegs. Since 'this' is not passed, we assume a global or context.
            // For translation purposes, we assume a global Parseador_t instance or similar context is available.
            // Here we mock the check.
            bool hasVolumes = false; // Mock
            
            if (hasVolumes) { 
                std::vector<ConformalMedia_t> conformal_media = buildConformalMedia(/*regs*/ {});
                std::vector<side_tris_map_t> side_to_triangles_maps = buildSideMaps(/*regs*/ {});
                for (m = 0; m < conformal_media.size(); m++) {
                    if (conformal_media[m].time_step_scale_factor < min_scale_factor) {
                        min_scale_factor = conformal_media[m].time_step_scale_factor;
                    }
                }
                dt = (1.0 / (cluz * std::sqrt(
                    (1.0 / sgg.DX) * (1.0 / sgg.DX) +
                    (1.0 / sgg.DY) * (1.0 / sgg.DY) +
                    (1.0 / sgg.DZ) * (1.0 / sgg.DZ)
                )));
                if (sgg.dt > dt * min_scale_factor) {
                    std::cout << "-- Conformal geometry requires a time step change" << std::endl;
                    std::cout << "Previous time step: " << sgg.dt << std::endl;
                    sgg.dt = dt * min_scale_factor;
                    std::cout << "New time step: " << sgg.dt << std::endl;
                }
            } else {
                std::vector<ConformalMedia_t> conformal_media(0);
            }
        }

        // Count media
        contamedia = 1;
        // Assuming 'this' has pmcregs. Mocking access.
        bool hasPMC = false; // Mock
        if (hasPMC) {
            contamedia = 2;
        }

        // Add other media counts (mocking 'this' members)
        contamedia += 10; // Placeholder for DielRegs, FRQDEPMATS, etc.
        contamedia += 7;  // LossyThinSurfs
        contamedia += 5;  // wires
        contamedia += 5;  // swires
        
        if (run_with_dmma) {
            contamedia += 10; // tSlots
        }

        if (MEDIOEXTRA.exists) {
            contamedia += 1;
            MEDIOEXTRA.index = contamedia;
        }

        contamedia += 2; // no_use, no_use_notouch
        contamedia += 1; // nodal sources

        std::vector<REAL> edge_ratios = getDifferentEdgeRatios({});
        std::vector<REAL> face_ratios = getDifferentFaceRatios({});
        contamedia += edge_ratios.size() + face_ratios.size();
        
        if (std::find(edge_ratios.begin(), edge_ratios.end(), 0.0) != edge_ratios.end()) {
            contamedia -= 1;
        }
        if (std::find(face_ratios.begin(), face_ratios.end(), 0.0) != face_ratios.end()) {
            contamedia -= 1;
        }

#ifdef CompileWithMTLN
        contamedia += 10; // mtln
#endif

        sgg.NumMedia = contamedia;
        sgg.AllocMed = contamedia;

        // Allocate media array
        sgg.Med.resize(contamedia + 1);

        // Set BoundingBox
        // Mocking sgg.Alloc access
        sgg.Alloc[iHx].XI = 0; sgg.Alloc[iHx].XE = 10;
        sgg.Alloc[iHy].YI = 0; sgg.Alloc[iHy].YE = 10;
        sgg.Alloc[iHz].ZI = 0; sgg.Alloc[iHz].ZE = 10;

        BoundingBox.XI = sgg.Alloc[iHx].XI;
        BoundingBox.XE = sgg.Alloc[iHx].XE;
        BoundingBox.YI = sgg.Alloc[iHy].YI;
        BoundingBox.YE = sgg.Alloc[iHy].YE;
        BoundingBox.ZI = sgg.Alloc[iHz].ZI;
        BoundingBox.ZE = sgg.Alloc[iHz].ZE;

        // Assign Alloc indices
        Alloc_iEx_XI = sgg.Alloc[iEx].XI; Alloc_iEx_XE = sgg.Alloc[iEx].XE;
        Alloc_iEx_YI = sgg.Alloc[iEx].YI; Alloc_iEx_YE = sgg.Alloc[iEx].YE;
        Alloc_iEx_ZI = sgg.Alloc[iEx].ZI; Alloc_iEx_ZE = sgg.Alloc[iEx].ZE;
        
        Alloc_iEy_XI = sgg.Alloc[iEy].XI; Alloc_iEy_XE = sgg.Alloc[iEy].XE;
        Alloc_iEy_YI = sgg.Alloc[iEy].YI; Alloc_iEy_YE = sgg.Alloc[iEy].YE;
        Alloc_iEy_ZI = sgg.Alloc[iEy].ZI; Alloc_iEy_ZE = sgg.Alloc[iEy].ZE;
        
        Alloc_iEz_XI = sgg.Alloc[iEz].XI; Alloc_iEz_XE = sgg.Alloc[iEz].XE;
        Alloc_iEz_YI = sgg.Alloc[iEz].YI; Alloc_iEz_YE = sgg.Alloc[iEz].YE;
        Alloc_iEz_ZI = sgg.Alloc[iEz].ZI; Alloc_iEz_ZE = sgg.Alloc[iEz].ZE;
        
        Alloc_iHx_XI = sgg.Alloc[iHx].XI; Alloc_iHx_XE = sgg.Alloc[iHx].XE;
        Alloc_iHx_YI = sgg.Alloc[iHx].YI; Alloc_iHx_YE = sgg.Alloc[iHx].YE;
        Alloc_iHx_ZI = sgg.Alloc[iHx].ZI; Alloc_iHx_ZE = sgg.Alloc[iHx].ZE;
        
        Alloc_iHy_XI = sgg.Alloc[iHy].XI; Alloc_iHy_XE = sgg.Alloc[iHy].XE;
        Alloc_iHy_YI = sgg.Alloc[iHy].YI; Alloc_iHy_YE = sgg.Alloc[iHy].YE;
        Alloc_iHy_ZI = sgg.Alloc[iHy].ZI; Alloc_iHy_ZE = sgg.Alloc[iHy].ZE;
        
        Alloc_iHz_XI = sgg.Alloc[iHz].XI; Alloc_iHz_XE = sgg.Alloc[iHz].XE;
        Alloc_iHz_YI = sgg.Alloc[iHz].YI; Alloc_iHz_YE = sgg.Alloc[iHz].YE;
        Alloc_iHz_ZI = sgg.Alloc[iHz].ZI; Alloc_iHz_ZE = sgg.Alloc[iHz].ZE;

        field = 1;
        numertag = 0;

        // Allocate media arrays
        // Helper lambda to allocate 3D vector
        auto alloc_3d = [&](std::vector<std::vector<std::vector<INT4>>>& vec, INT4 x1, INT4 x2, INT4 y1, INT4 y2, INT4 z1, INT4 z2) {
            vec.resize(x2 - x1 + 1);
            for (INT4 x = 0; x < vec.size(); ++x) {
                vec[x].resize(y2 - y1 + 1);
                for (INT4 y = 0; y < vec[x].size(); ++y) {
                    vec[x][y].resize(z2 - z1 + 1, 0);
                }
            }
        };

        alloc_3d(media.sggMtag, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHz_ZI, Alloc_iHz_ZE);
        alloc_3d(media.sggMiNo, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHz_ZI, Alloc_iHz_ZE);

        alloc_3d(tag_numbers.edge.x, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE);
        alloc_3d(tag_numbers.edge.y, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE);
        alloc_3d(tag_numbers.edge.z, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE);

        alloc_3d(tag_numbers.face.x, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE);
        alloc_3d(tag_numbers.face.y, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE);
        alloc_3d(tag_numbers.face.z, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE);

        alloc_3d(media.sggMiEx, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE);
        alloc_3d(media.sggMiEy, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE);
        alloc_3d(media.sggMiEz, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE);
        alloc_3d(media.sggMiHx, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE);
        alloc_3d(media.sggMiHy, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE);
        alloc_3d(media.sggMiHz, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE);

        // Initialize arrays to 0 or 1
        auto fill_3d = [&](std::vector<std::vector<std::vector<INT4>>>& vec, INT4 val) {
            for (auto& x : vec) {
                for (auto& y : x) {
                    for (auto& z : y) {
                        z = val;
                    }
                }
            }
        };

        fill_3d(media.sggMtag, 0);
        fill_3d(tag_numbers.edge.x, 0);
        fill_3d(tag_numbers.edge.y, 0);
        fill_3d(tag_numbers.edge.z, 0);
        fill_3d(tag_numbers.face.x, 0);
        fill_3d(tag_numbers.face.y, 0);
        fill_3d(tag_numbers.face.z, 0);

        fill_3d(media.sggMiNo, 1);
        fill_3d(media.sggMiEx, 1);
        fill_3d(media.sggMiEy, 1);
        fill_3d(media.sggMiEz, 1);
        fill_3d(media.sggMiHx, 1);
        fill_3d(media.sggMiHy, 1);
        fill_3d(media.sggMiHz, 1);

        // Plane Waves
        INT4 tama = 1; // Mock this->plnSrc.nc
        amplitud = 1.0;
        sgg.NumPlaneWaves = tama;
        sgg.PlaneWave.resize(tama);

        for (i = 1; i <= sgg.NumPlaneWaves; ++i) {
            // Mock coordinates
            punto.XI = 0; punto.XE = 10;
            punto.YI = 0; punto.YE = 10;
            punto.ZI = 0; punto.ZE = 10;

            // Adjustments
            if (punto.XI == SINPML_fullsize[iHx].bounds.XI) punto.XI -= 5;
            if (punto.XE == SINPML_fullsize[iHx].bounds.XE) punto.XE += 5;
            if (punto.YI == SINPML_fullsize[iHy].bounds.YI) punto.YI -= 5;
            if (punto.YE == SINPML_fullsize[iHy].bounds.YE) punto.YE += 5;
            if (punto.ZI == SINPML_fullsize[iHz].bounds.ZI) punto.ZI -= 5;
            if (punto.ZE == SINPML_fullsize[iHz].bounds.ZE) punto.ZE += 5;

            sgg.PlaneWave[i-1].isRC = false; // Mock
            sgg.PlaneWave[i-1].numModes = 1; // Mock
            sgg.PlaneWave[i-1].incertmax = 0.0; // Mock

            if (sgg.PlaneWave[i-1].isRC) {
                INT4 nm = sgg.PlaneWave[i-1].numModes;
                sgg.PlaneWave[i-1].px.resize(nm, 0.0);
                sgg.PlaneWave[i-1].py.resize(nm, 0.0);
                sgg.PlaneWave[i-1].pz.resize(nm, 0.0);
                sgg.PlaneWave[i-1].ex.resize(nm, 0.0);
                sgg.PlaneWave[i-1].ey.resize(nm, 0.0);
                sgg.PlaneWave[i-1].ez.resize(nm, 0.0);
                sgg.PlaneWave[i-1].INCERT.resize(nm, 0.0);
                
                dummy_px.resize(nm);
                dummy_py.resize(nm);
                dummy_pz.resize(nm);
                dummy_ex.resize(nm);
                dummy_ey.resize(nm);
                dummy_ez.resize(nm);
                dummy_INCERT.resize(nm);

                if (layoutnumber == 0) {
                    populatePlaneWaveRC(sgg.PlaneWave[i-1], simu_devia);
                }

#ifdef CompileWithMPI
                MPI_Barrier(MPI_COMM_WORLD, &IERR);
                MPI_AllReduce(sgg.PlaneWave[i-1].px.data(), dummy_px.data(), nm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &IERR);
                MPI_AllReduce(sgg.PlaneWave[i-1].py.data(), dummy_py.data(), nm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &IERR);
                MPI_AllReduce(sgg.PlaneWave[i-1].pz.data(), dummy_pz.data(), nm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &IERR);
                MPI_AllReduce(sgg.PlaneWave[i-1].ex.data(), dummy_ex.data(), nm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &IERR);
                MPI_AllReduce(sgg.PlaneWave[i-1].ey.data(), dummy_ey.data(), nm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &IERR);
                MPI_AllReduce(sgg.PlaneWave[i-1].ez.data(), dummy_ez.data(), nm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &IERR);
                MPI_AllReduce(sgg.PlaneWave[i-1].INCERT.data(), dummy_INCERT.data(), nm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &IERR);
                MPI_Barrier(MPI_COMM_WORLD, &IERR);

                sgg.PlaneWave[i-1].px = dummy_px;
                sgg.PlaneWave[i-1].py = dummy_py;
                sgg.PlaneWave[i-1].pz = dummy_pz;
                sgg.PlaneWave[i-1].ex = dummy_ex;
                sgg.PlaneWave[i-1].ey = dummy_ey;
                sgg.PlaneWave[i-1].ez = dummy_ez;
                sgg.PlaneWave[i-1].INCERT = dummy_INCERT;
#endif
                dummy_px.clear();
                dummy_py.clear();
                dummy_pz.clear();
                dummy_ex.clear();
                dummy_ey.clear();
                dummy_ez.clear();
                dummy_INCERT.clear();
            } else {
                sgg.PlaneWave[i-1].numModes = 1;
                sgg.PlaneWave[i-1].px.resize(1, 0.0);
                sgg.PlaneWave[i-1].py.resize(1, 0.0);
                sgg.PlaneWave[i-1].pz.resize(1, 0.0);
                sgg.PlaneWave[i-1].ex.resize(1, 0.0);
                sgg.PlaneWave[i-1].ey.resize(1, 0.0);
                sgg.PlaneWave[i-1].ez.resize(1, 0.0);
                sgg.PlaneWave[i-1].INCERT.resize(1, 0.0);
            }
        }
    }

    void read_limits_nogeom() {
        /* TODO: Read geometry limits without full geometry parsing. From preprocess_geom.F90 */
    }

    void AssigLossyOrPECtoNodes() {
        /* TODO: Assign lossy or PEC material to nodes. From preprocess_geom.F90 */
    }

    void searchtag() {
        /* TODO: Search for tag definitions. From preprocess_geom.F90 */
    }

    void checkDielectricTagForDuplicate() {
        /* TODO: Check for duplicate dielectric tags. From preprocess_geom.F90 */
    }

    void checkAnimatedTagForDuplicate() {
        /* TODO: Check for duplicate animated tags. From preprocess_geom.F90 */
    }

    void checkLossyTagForDuplicate() {
        /* TODO: Check for duplicate lossy tags. From preprocess_geom.F90 */
    }

} // namespace Preprocess_m

allocate(sgg%PlaneWave(i)%ex(1:sgg%PlaneWave(i)%numModes))
            allocate(sgg%PlaneWave(i)%ey(1:sgg%PlaneWave(i)%numModes))
            allocate(sgg%PlaneWave(i)%ez(1:sgg%PlaneWave(i)%numModes))
            allocate(sgg%PlaneWave(i)%INCERT(1:sgg%PlaneWave(i)%numModes))
            !
            ez = amplitud * Cos(this%plnSrc%collection(i)%alpha)
            ey = amplitud * Sin(this%plnSrc%collection(i)%alpha) * Sin(this%plnSrc%collection(i)%beta)
            ex = amplitud * Sin(this%plnSrc%collection(i)%alpha) * Cos(this%plnSrc%collection(i)%beta)
            pz = Cos(this%plnSrc%collection(i)%theta)
            py = Sin(this%plnSrc%collection(i)%theta) * Sin(this%plnSrc%collection(i)%phi)
            px = Sin(this%plnSrc%collection(i)%theta) * Cos(this%plnSrc%collection(i)%phi)
            !ojo con estos redondeos.
            !!!if (Abs(ex/amplitud) < 1e-4) ex = 0.0_RKIND
            !!!if (Abs(ey/amplitud) < 1e-4) ey = 0.0_RKIND
            !!!if (Abs(ez/amplitud) < 1e-4) ez = 0.0_RKIND
            !!!if (Abs(px) < 1e-4) px = 0.0_RKIND
            !!!if (Abs(py) < 1e-4) py = 0.0_RKIND
            !!!if (Abs(pz) < 1e-4) pz = 0.0_RKIND
            if (Abs(px*ex+py*ey+pz*ez) >= 1e-4) then
               write(buff,*) 'NO TEM PLANEWAVE',ex,ey,ez,px,py,pz,(px*ex+py*ey+pz*ez),this%plnSrc%collection(i)%alpha,  &
                  this%plnSrc%collection(i)%beta,this%plnSrc%collection(i)%theta,this%plnSrc%collection(i)%phi
               call STOPONERROR(layoutnumber,num_procs,buff)
            end if
            !
            sgg%PlaneWave(i)%px(1) = px
            sgg%PlaneWave(i)%py(1) = py
            sgg%PlaneWave(i)%pz(1) = pz
            sgg%PlaneWave(i)%ex(1) = ex
            sgg%PlaneWave(i)%ey(1) = ey
            sgg%PlaneWave(i)%ez(1) = ez
            sgg%PlaneWave(i)%INCERT(1)=0.0_RKIND
         end if
         sgg%PlaneWave(i)%fichero%name = trim(adjustl(this%plnSrc%collection(i)%nombre_fichero))
         sgg%PlaneWave(i)%esqx1 = Min(punto%XI, punto%XE)
         sgg%PlaneWave(i)%esqy1 = Min(punto%YI, punto%YE)
         sgg%PlaneWave(i)%esqz1 = Min(punto%ZI, punto%ZE)
         sgg%PlaneWave(i)%esqx2 = Max(punto%XI, punto%XE)
         sgg%PlaneWave(i)%esqy2 = Max(punto%YI, punto%YE)
         sgg%PlaneWave(i)%esqz2 = Max(punto%ZI, punto%ZE)
      end do
      !Media parsing
      !Default
      !background
      sgg%Med%Priority = prior_BV
      sgg%Med%Epr = 1.0
      sgg%Med%Sigma = 0.0
      sgg%Med%Sigmareasignado = .false. !solo afecta a un chequeo de errores en lumped 120123
      sgg%Med%Mur = 1.0
      sgg%Med%SigmaM = 0.0
      sgg%Med%Is%Interfase = .FALSE.
      sgg%Med%Is%PMLbody = .false.
      sgg%Med%Is%Needed = .TRUE.
      sgg%Med%Is%Anisotropic = .FALSE.
      sgg%Med%Is%Dielectric = .FALSE.
      sgg%Med%Is%EDispersive = .FALSE.
      sgg%Med%Is%EDispersiveAnis = .FALSE.
      sgg%Med%Is%MDispersive = .FALSE.
      sgg%Med%Is%MDispersiveAnis = .FALSE.
      sgg%Med%Is%Lumped = .FALSE.
      sgg%Med%Is%SGBC = .FALSE.
      sgg%Med%Is%SGBCDispersive = .FALSE.
      sgg%Med%Is%Lossy = .FALSE.
      sgg%Med%Is%multiport = .FALSE.
      sgg%Med%Is%multiportpadding = .FALSE.
      sgg%Med%Is%AnisMultiport = .FALSE.
      sgg%Med%Is%ThinWire = .FALSE.
      sgg%Med%Is%Multiwire = .FALSE.
      sgg%Med%Is%SlantedWire = .FALSE.
      sgg%Med%Is%ThinSlot = .FALSE.
      sgg%Med%Is%PEC = .FALSE.
      sgg%Med%Is%ConformalPEC = .FALSE.
      sgg%Med%Is%PMC = .FALSE.
      sgg%Med%Is%PML = .FALSE.
      sgg%Med%Is%Volume = .FALSE.
      sgg%Med%Is%Surface = .FALSE.
      sgg%Med%Is%Line = .FALSE.
      sgg%Med%Is%already_YEEadvanced_byconformal = .FALSE.
      sgg%Med%Is%split_and_useless = .FALSE.
      !ojo tocar tambien en el readjust de healing si se crean nuevos flags
      !
      !medio PEC y PML es intrascendente si es surface o volume
      !son los de prioridad mas alta y siempre contienen a sus campos tangenciales electricos
      !Background    only differences from default are needed
      sgg%Med(1)%Priority = prior_BV
      sgg%Med(1)%Epr = this%mats%mats(1)%eps / Eps0
      sgg%Med(1)%Sigma = this%mats%mats(1)%Sigma
      sgg%Med(1)%Mur = this%mats%mats(1)%mu / Mu0
      sgg%Med(1)%SigmaM = this%mats%mats(1)%SigmaM
      sgg%Med(1)%Is%Dielectric = .false. !considero el vacio como NO dielectrico '251114
      sgg%Med(1)%Is%Volume = .false.  !considero el vacio como no volumic false '251114
      !
      sgg%Med(0)%Is%PEC = .TRUE.
      sgg%Med(0)%Is%Needed = .TRUE.
      sgg%Med(0)%Priority = prior_PEC
      sgg%Med(0)%Epr = this%mats%mats(1)%eps / Eps0
      sgg%Med(0)%Sigma = 1.0e29_RKIND
      sgg%Med(0)%Mur = this%mats%mats(1)%mu / Mu0
      sgg%Med(0)%SigmaM = 0.0_RKIND
      !CAPA EXTRA
      !Background    only differences from default are needed

      if (medioextra%exists) then
         !!!!estimate in terms of percentage of the maximum PML conductivity the conductivity of the extra medium
         !!!This info is available from read_limits_nogeom
         !the calculus is taken from borderscpml.F90
         sig_max=0.0_RKIND
         do o=1,3
            do p=1,2
               if ((o == 1).and.(p == 1)) del=sgg%dx(SINPML_fullsize(iHx)%XI)
               if ((o == 1).and.(p == 2)) del=sgg%dx(SINPML_fullsize(iHx)%XE-1)
               if ((o == 2).and.(p == 1)) del=sgg%dy(SINPML_fullsize(iHy)%YI)
               if ((o == 2).and.(p == 2)) del=sgg%dy(SINPML_fullsize(iHy)%YE-1)
               if ((o == 3).and.(p == 1)) del=sgg%dz(SINPML_fullsize(iHz)%ZI)
               if ((o == 3).and.(p == 2)) del=sgg%dz(SINPML_fullsize(iHz)%ZE-1)
               if (sgg%PML%NumLayers(o,p) /= 0) then
                  if ((sgg%PML%NumLayers(o,p) == 10).or.(sgg%PML%NumLayers(o,p) == 5)) then
                     sig_max = max(sig_max, 0.8*(sgg%PML%orden(o,p)+1)/(zvac*del))
                  else
                     if (sgg%PML%CoeffReflPML(o,p)==1.0_RKIND) then
                        !realmente en el borderscpml
                        !sig_max(sig_max,-((log( 0.99999d0                 )*(sgg%PML%orden(o,p)+1))/ &
                        !    (2.0_RKIND *sqrt(Mu0/eps0)*sgg%PML%NumLayers(o,p)*del)))
                        !trampa para que entonces tome la conductividad autentica que se especifique y poder anular las PML y solo dejar capa fisica !!?!?
                        sig_max = 1.0_RKIND
                     else
                        sig_max = max(sig_max,-((log(sgg%PML%CoeffReflPML(o,p))*(sgg%PML%orden(o,p)+1))/ &
                           (2.0_RKIND *sqrt(Mu0/eps0)*sgg%PML%NumLayers(o,p)*del)))
                     end if
                  end if
               end if
            end do
         end do
         MEDIOEXTRA%sigma = MEDIOEXTRA%sigma * sig_max !la especificacion se da en terminos de tanto por uno en la linea de comandos
         !
         sgg%Med(MEDIOEXTRA%index)%Epr = this%mats%mats(1)%eps / Eps0 !luego se machaca este valor
         sgg%Med(MEDIOEXTRA%index)%Sigma = MEDIOEXTRA%sigma !luego se machaca este valor
         sgg%Med(MEDIOEXTRA%index)%Mur = this%mats%mats(1)%mu / Mu0 !luego se machaca este valor
         sgg%Med(MEDIOEXTRA%index)%SigmaM = 0.0_RKIND !solo lo creo para las tangenciales electricas
         sgg%Med(MEDIOEXTRA%index)%Priority = prior_PEC
         sgg%Med(MEDIOEXTRA%index)%Is%Dielectric = .TRUE.
         sgg%Med(MEDIOEXTRA%index)%Is%Volume = .TRUE.
         sgg%Med(MEDIOEXTRA%index)%Is%PML = .TRUE.
      end if
      !
      !barre los medios
      !Primero todos los pec
      !PECRegions
      !volumenes
      !el medio 0 se reserva para PEC
      !regiones PEC
      !
      if ((this%pecregs%nvols)+(this%pecregs%nsurfs)+(this%pecregs%nLINS) /= 0) then
         pecmedio = 0
         tama = (this%pecregs%nvols)
         !BODYes
         do i = 1, tama
            punto%XI = this%pecregs%vols(i)%XI
            punto%XE = this%pecregs%vols(i)%XE
            punto%YI = this%pecregs%vols(i)%YI
            punto%YE = this%pecregs%vols(i)%YE
            punto%ZI = this%pecregs%vols(i)%ZI
            punto%ZE = this%pecregs%vols(i)%ZE
            numertag = searchtag(tagtype,this%pecregs%vols(i)%tag)
            call CreateVolumeMM(layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, pecmedio)
         end do
         !SURFs
         tama = (this%pecregs%nsurfs)
         do i = 1, tama
            punto%XI = this%pecregs%surfs(i)%XI
            punto%XE = this%pecregs%surfs(i)%XE
            punto%YI = this%pecregs%surfs(i)%YI
            punto%YE = this%pecregs%surfs(i)%YE
            punto%ZI = this%pecregs%surfs(i)%ZI
            punto%ZE = this%pecregs%surfs(i)%ZE
            orientacion = this%pecregs%surfs(i)%or
            numertag = searchtag(tagtype,this%pecregs%surfs(i)%tag)
            call CreateSurfaceMM(layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, &
            & Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, &
            & Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
            & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, &
            & Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, &
            & Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, &
            & Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, &
            & sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, pecmedio)
         end do
         !LINs
         tama = (this%pecregs%nLINS)
         do i = 1, tama
            punto%XI = this%pecregs%lins(i)%XI
            punto%XE = this%pecregs%lins(i)%XE
            punto%YI = this%pecregs%lins(i)%YI
            punto%YE = this%pecregs%lins(i)%YE
            punto%ZI = this%pecregs%lins(i)%ZI
            punto%ZE = this%pecregs%lins(i)%ZE
            orientacion = this%pecregs%lins(i)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%pecregs%lins(i)%tag)
            call CreateLineMM(layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & pecmedio, isathinwire, verbose,numeroasignaciones)
         end do
         !regiones PEC
      end if

      !el medio 1 se reserva para sustrato  y saltamos
      contamedia = 1

      !para el conformal !debe ser tipicamente contamedia =1+1=2 pq el 0 es pec y el 1 es vacio. Ojo cambiado de sitio el PMC porque podia hacer que fuesen 3 y 4. 130220!!! y puede haber error pq por ahi se comprueba el 2 y el 3
      contamedia = contamedia + 1
      sgg%Med(contamedia)%Is%already_YEEadvanced_byconformal = .TRUE.
      !debe ser contamedia =2+1=3
      contamedia = contamedia + 1
      sgg%Med(contamedia)%Is%split_and_useless = .TRUE.

!!!!cambiado aqui 130220

      !materialList
      !regiones PMC
      if ((this%pmcregs%nvols)+(this%pmcregs%nsurfs)+(this%pmcregs%nLINS) /= 0) then
         !los PMC de existir tienen todos indice 2
         contamedia =contamedia+1      !!!!contamedia = 2 !!!ufff. cambiado a 130220 por posible bug con conformal si algun dia habia regiones PMC
         sgg%Med(contamedia)%Epr = sgg%Med(1)%Epr
         sgg%Med(contamedia)%Mur = sgg%Med(1)%Mur
         sgg%Med(contamedia)%Sigma = 0.0_RKIND
         sgg%Med(contamedia)%SigmaM = 1.0e29_RKIND
         sgg%Med(contamedia)%Priority = prior_PMC
         sgg%Med(contamedia)%Is%PMC = .TRUE.
         !BODYes
         tama = (this%pmcregs%nvols)
         do i = 1, tama
            punto%XI = this%pmcregs%vols(i)%XI
            punto%XE = this%pmcregs%vols(i)%XE
            punto%YI = this%pmcregs%vols(i)%YI
            punto%YE = this%pmcregs%vols(i)%YE
            punto%ZI = this%pmcregs%vols(i)%ZI
            punto%ZE = this%pmcregs%vols(i)%ZE
            !
            !
            numertag = searchtag(tagtype,this%pmcregs%vols(i)%tag)
            call CreateVolumeMM(layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         end do
         !SURFs
         tama = (this%pmcregs%nsurfs)
         do i = 1, tama
            punto%XI = this%pmcregs%surfs(i)%XI
            punto%XE = this%pmcregs%surfs(i)%XE
            punto%YI = this%pmcregs%surfs(i)%YI
            punto%YE = this%pmcregs%surfs(i)%YE
            punto%ZI = this%pmcregs%surfs(i)%ZI
            punto%ZE = this%pmcregs%surfs(i)%ZE
            orientacion = this%pmcregs%surfs(i)%or
            numertag = searchtag(tagtype,this%pmcregs%surfs(i)%tag)
            call CreateSurfaceMM(layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia)
         end do
         !LINs
         tama = (this%pmcregs%nLINS)
         do i = 1, tama
            punto%XI = this%pmcregs%lins(i)%XI
            punto%XE = this%pmcregs%lins(i)%XE
            punto%YI = this%pmcregs%lins(i)%YI
            punto%YE = this%pmcregs%lins(i)%YE
            punto%ZI = this%pmcregs%lins(i)%ZI
            punto%ZE = this%pmcregs%lins(i)%ZE
            orientacion = this%pmcregs%lins(i)%or
            isathinwire = .FALSE.
            numertag = searchtag(tagtype,this%pmcregs%lins(i)%tag)
            call CreateLineMM(layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire,verbose,numeroasignaciones)
            !
         end do
         !fin regions PMC
      end if
!!!!fin cambiado 130220

      !NonMetalREgions
      !BODYes
      tama = (this%DielRegs%nvols)
      do i = 1, tama
         contamedia = contamedia + 1
         sgg%Med(contamedia)%Is%Dielectric = .TRUE.
         sgg%Med(contamedia)%Priority = prior_IB
         sgg%Med(contamedia)%Epr = this%DielRegs%vols(i)%eps / Eps0
         sgg%Med(contamedia)%Sigma = this%DielRegs%vols(i)%Sigma
         sgg%Med(contamedia)%Mur =    this%DielRegs%vols(i)%mu / Mu0
         sgg%Med(contamedia)%SigmaM = this%DielRegs%vols(i)%SigmaM
!!!!pmlbody
         if (this%DielRegs%vols(i)%PMLbody) then
            sgg%Med(contamedia)%Priority = prior_pmlbody !machaca con una prioridad superior a la de thin wires y backgroud !prueba HOLD 251019 coax
            sgg%Med(contamedia)%Is%PMLbody = .true.
           allocate(sgg%Med(contamedia)%PMLbody(1))
            sgg%Med(contamedia)%PMLbody(1)%orient    = this%DielRegs%vols(i)%orient
         end if
!!!!!
         tama2 = (this%DielRegs%vols(i)%n_c2P)
         do j = 1, tama2
            if ((J==1).and.(this%DielRegs%vols(i)%PMLbody)) sgg%Med(contamedia)%PMLbody(1)%orient = this%DielRegs%vols(i)%c2P(j)%OR !ES IGUAL PARA TODOS
            punto%XI = this%DielRegs%vols(i)%c2P(j)%XI
            punto%XE = this%DielRegs%vols(i)%c2P(j)%XE
            punto%YI = this%DielRegs%vols(i)%c2P(j)%YI
            punto%YE = this%DielRegs%vols(i)%c2P(j)%YE
            punto%ZI = this%DielRegs%vols(i)%c2P(j)%ZI
            punto%ZE = this%DielRegs%vols(i)%c2P(j)%ZE
            numertag = searchtag(tagtype,this%DielRegs%vols(i)%c2P(j)%tag)
            call CreateVolumeMM(layoutnumber, media%sggMtag, tag_numbers, numertag, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
            & media%sggMiHx, media%sggMiHy, media%sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg%Med, sgg%NumMedia, sgg%EShared, BoundingBox, punto, contamedia)
         end do
         tama3 = (this%DielRegs%vols(i)%n_c1P)
         do j = 1, tama3
            if ((J==1).and.(this%DielRegs%vols(i)%PMLbody)) sgg%Med(contamedia)%PMLbody(1)%orient = this%DielRegs%vols(i)%c1P(j)%OR !ES IGUAL PARA TODOS
            punto%XI = this%DielRegs%vols(i)%c1P(j)%XI
            punto%XE = this%DielRegs%vols(i)%c1P(j)%XI
            punto%YI = this%DielRegs%vols(i)%c1P(j)%YI
            punto%YE = this%DielRegs%vols(i)%c1P(j)%YI

punto.ZI = this->DielRegs.vols(i).c1P(j).ZI;
            punto.ZE = this->DielRegs.vols(i).c1P(j).ZI;
            
            numertag = searchtag(tagtype, this->DielRegs.vols(i).c1P(j).tag);
            CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                           media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                           Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                           Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                           Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                           Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                           Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
         }
      }
      // SURFs
      tama = (this->DielRegs.nsurfs);
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         sgg.Med(contamedia).Is.Dielectric = true;
         sgg.Med(contamedia).Priority = prior_IS;
         sgg.Med(contamedia).Epr = this->DielRegs.surfs(i).eps / Eps0;
         sgg.Med(contamedia).Sigma = this->DielRegs.surfs(i).Sigma;
         sgg.Med(contamedia).Mur = this->DielRegs.surfs(i).mu / Mu0;
         sgg.Med(contamedia).SigmaM = this->DielRegs.surfs(i).SigmaM;
         tama2 = (this->DielRegs.surfs(i).n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->DielRegs.surfs(i).c2P(j).XI;
            punto.XE = this->DielRegs.surfs(i).c2P(j).XE;
            punto.YI = this->DielRegs.surfs(i).c2P(j).YI;
            punto.YE = this->DielRegs.surfs(i).c2P(j).YE;
            punto.ZI = this->DielRegs.surfs(i).c2P(j).ZI;
            punto.ZE = this->DielRegs.surfs(i).c2P(j).ZE;
            orientacion = this->DielRegs.surfs(i).c2P(j).or;
            numertag = searchtag(tagtype, this->DielRegs.surfs(i).c2P(j).tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
         }
         tama3 = (this->DielRegs.surfs(i).n_c1P);

         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->DielRegs.surfs(i).c1P(j).XI;
            punto.XE = this->DielRegs.surfs(i).c1P(j).XI;
            punto.YI = this->DielRegs.surfs(i).c1P(j).YI;
            punto.YE = this->DielRegs.surfs(i).c1P(j).YI;
            punto.ZI = this->DielRegs.surfs(i).c1P(j).ZI;
            punto.ZE = this->DielRegs.surfs(i).c1P(j).ZI;
            orientacion = this->DielRegs.surfs(i).c1P(j).or;
            numertag = searchtag(tagtype, this->DielRegs.surfs(i).c1P(j).tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
         }
      }
      // LINs
      tama = (this->DielRegs.nLINS);
      for (i = 1; i <= tama; ++i) {
         numeroasignaciones = 0; // solo lo usa lumped para echarselo al primer y el resto ponerlo a PEC
         contamedia = contamedia + 1;
         sgg.Med(contamedia).Is.Dielectric = true;
         sgg.Med(contamedia).Priority = prior_IL;
         sgg.Med(contamedia).Epr = this->DielRegs.lins(i).eps / Eps0;
         sgg.Med(contamedia).Sigma = this->DielRegs.lins(i).Sigma;
         sgg.Med(contamedia).Mur = this->DielRegs.lins(i).mu / Mu0;
         sgg.Med(contamedia).SigmaM = this->DielRegs.lins(i).SigmaM;
         // lumped
         if (this->DielRegs.lins(i).resistor) {
            sgg.Med(contamedia).Is.Lumped = true;
            sgg.Med(contamedia).Is.lossy = true; // importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med(contamedia).lumped.resize(1);
            sgg.Med(contamedia).lumped[0].resistor = true;
            sgg.Med(contamedia).lumped[0].inductor = false;
            sgg.Med(contamedia).lumped[0].capacitor = false;
            sgg.Med(contamedia).lumped[0].diodo = false;
            sgg.Med(contamedia).lumped[0].R = this->DielRegs.lins(i).R;
            sgg.Med(contamedia).lumped[0].L = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].C = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].R_devia = this->DielRegs.lins(i).R_devia;
            sgg.Med(contamedia).lumped[0].L_devia = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].C_devia = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].Rtime_on = this->DielRegs.lins(i).Rtime_on;
            sgg.Med(contamedia).lumped[0].Rtime_off = this->DielRegs.lins(i).Rtime_off;
            sgg.Med(contamedia).lumped[0].DiodB = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].DiodIsat = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].orient = this->DielRegs.lins(i).DiodOri;
         } else if (this->DielRegs.lins(i).inductor) {
            sgg.Med(contamedia).Is.Lumped = true;
            sgg.Med(contamedia).Is.lossy = true; // importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med(contamedia).lumped.resize(1);
            sgg.Med(contamedia).lumped[0].resistor = false;
            sgg.Med(contamedia).lumped[0].inductor = true;
            sgg.Med(contamedia).lumped[0].capacitor = false;
            sgg.Med(contamedia).lumped[0].diodo = false;
            sgg.Med(contamedia).lumped[0].R = this->DielRegs.lins(i).R;
            sgg.Med(contamedia).lumped[0].L = this->DielRegs.lins(i).L;
            sgg.Med(contamedia).lumped[0].C = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].R_devia = this->DielRegs.lins(i).R_devia;
            sgg.Med(contamedia).lumped[0].L_devia = this->DielRegs.lins(i).L_devia;
            sgg.Med(contamedia).lumped[0].C_devia = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].Rtime_on = 0.0; // irrelevant
            sgg.Med(contamedia).lumped[0].Rtime_off = 0.0; // irrelevant
            sgg.Med(contamedia).lumped[0].DiodB = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].DiodIsat = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].orient = this->DielRegs.lins(i).DiodOri;
         } else if (this->DielRegs.lins(i).capacitor) {
            sgg.Med(contamedia).Is.Lumped = true;
            sgg.Med(contamedia).Is.lossy = true; // importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med(contamedia).lumped.resize(1);
            sgg.Med(contamedia).lumped[0].resistor = false;
            sgg.Med(contamedia).lumped[0].inductor = false;
            sgg.Med(contamedia).lumped[0].capacitor = true;
            sgg.Med(contamedia).lumped[0].diodo = false;
            sgg.Med(contamedia).lumped[0].R = this->DielRegs.lins(i).R;
            sgg.Med(contamedia).lumped[0].L = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].C = this->DielRegs.lins(i).C;
            sgg.Med(contamedia).lumped[0].R_devia = this->DielRegs.lins(i).R_devia;
            sgg.Med(contamedia).lumped[0].L_devia = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].C_devia = this->DielRegs.lins(i).C_devia;
            sgg.Med(contamedia).lumped[0].Rtime_on = 0.0; // irrelevant
            sgg.Med(contamedia).lumped[0].Rtime_off = 0.0; // irrelevant
            sgg.Med(contamedia).lumped[0].DiodB = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].DiodIsat = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].orient = this->DielRegs.lins(i).DiodOri;
         } else if (this->DielRegs.lins(i).diodo) {
            // 27/08/15 diodos aun no soportados
            std::string buff = "Lumped Diodes currently unsupported. .";
            STOPONERROR(layoutnumber, num_procs, buff);
            
            sgg.Med(contamedia).Is.Lumped = true;
            sgg.Med(contamedia).Is.lossy = true; // importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med(contamedia).lumped.resize(1);
            sgg.Med(contamedia).lumped[0].resistor = false;
            sgg.Med(contamedia).lumped[0].inductor = false;
            sgg.Med(contamedia).lumped[0].capacitor = false;
            sgg.Med(contamedia).lumped[0].diodo = true;
            sgg.Med(contamedia).lumped[0].R = this->DielRegs.lins(i).R;
            sgg.Med(contamedia).lumped[0].Rtime_on = 0.0; // irrelevant
            sgg.Med(contamedia).lumped[0].Rtime_off = 0.0; // irrelevant
            sgg.Med(contamedia).lumped[0].L = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].C = 0.0_RKIND;
            sgg.Med(contamedia).lumped[0].DiodB = this->DielRegs.lins(i).DiodB;
            sgg.Med(contamedia).lumped[0].DiodIsat = this->DielRegs.lins(i).DiodIsat;
            sgg.Med(contamedia).lumped[0].orient = this->DielRegs.lins(i).DiodOri;
         } else {
            sgg.Med(contamedia).Is.Lumped = false;
            if (!this->DielRegs.lins(i).plain) {
               std::string buff = "Buggy error 1 in preprocess lumped. .";
               STOPONERROR(layoutnumber, num_procs, buff);
            }
         }
         // fin lumped
         tama2 = (this->DielRegs.lins(i).n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->DielRegs.lins(i).c2P(j).XI;
            punto.XE = this->DielRegs.lins(i).c2P(j).XE;
            punto.YI = this->DielRegs.lins(i).c2P(j).YI;
            punto.YE = this->DielRegs.lins(i).c2P(j).YE;
            punto.ZI = this->DielRegs.lins(i).c2P(j).ZI;
            punto.ZE = this->DielRegs.lins(i).c2P(j).ZE;
            orientacion = this->DielRegs.lins(i).c2P(j).or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->DielRegs.lins(i).c2P(j).tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
         }
         tama3 = (this->DielRegs.lins(i).n_c1P);
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->DielRegs.lins(i).c1P(j).XI;
            punto.XE = this->DielRegs.lins(i).c1P(j).XI;
            punto.YI = this->DielRegs.lins(i).c1P(j).YI;
            punto.YE = this->DielRegs.lins(i).c1P(j).YI;
            punto.ZI = this->DielRegs.lins(i).c1P(j).ZI;
            punto.ZE = this->DielRegs.lins(i).c1P(j).ZI;
            orientacion = this->DielRegs.lins(i).c1P(j).or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->DielRegs.lins(i).c1P(j).tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
         }
      }


      // Anisotropic materials
      // materialList
      // BODYes
      tama = (this->ANIMATS.nvols);
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         sgg.Med(contamedia).Anisotropic.resize(1);
         sgg.Med(contamedia).Is.Anisotropic = true;
         sgg.Med(contamedia).Priority = prior_AB;
         sgg.Med(contamedia).Anisotropic[0].Epr = this->ANIMATS.vols(i).eps / Eps0;
         sgg.Med(contamedia).Anisotropic[0].Sigma = this->ANIMATS.vols(i).Sigma;
         sgg.Med(contamedia).Anisotropic[0].Mur = this->ANIMATS.vols(i).mu / Mu0;
         sgg.Med(contamedia).Anisotropic[0].SigmaM = this->ANIMATS.vols(i).SigmaM;
         tama2 = (this->ANIMATS.vols(i).n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->ANIMATS.vols(i).c2P(j).XI;
            punto.XE = this->ANIMATS.vols(i).c2P(j).XE;
            punto.YI = this->ANIMATS.vols(i).c2P(j).YI;
            punto.YE = this->ANIMATS.vols(i).c2P(j).YE;
            punto.ZI = this->ANIMATS.vols(i).c2P(j).ZI;
            punto.ZE = this->ANIMATS.vols(i).c2P(j).ZE;
            numertag = searchtag(tagtype, this->ANIMATS.vols(i).c2P(j).tag);
            CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                           media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                           Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                           Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                           Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                           Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                           Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
         }
         tama3 = (this->ANIMATS.vols(i).n_c1P);
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->ANIMATS.vols(i).c1P(j).XI;
            punto.XE = this->ANIMATS.vols(i).c1P(j).XI;
            punto.YI = this->ANIMATS.vols(i).c1P(j).YI;
            punto.YE = this->ANIMATS.vols(i).c1P(j).YI;
            punto.ZI = this->ANIMATS.vols(i).c1P(j).ZI;
            punto.ZE = this->ANIMATS.vols(i).c1P(j).ZI;
            
            numertag = searchtag(tagtype, this->ANIMATS.vols(i).c1P(j).tag);
            CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                           media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                           Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                           Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                           Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                           Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                           Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
         }
      }
      // SURFs
      tama = (this->ANIMATS.nsurfs);
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         sgg.Med(contamedia).Anisotropic.resize(1);
         sgg.Med(contamedia).Is.Anisotropic = true;
         sgg.Med(contamedia).Priority = prior_IS;
         sgg.Med(contamedia).Anisotropic[0].Epr = this->ANIMATS.surfs(i).eps / Eps0;
         sgg.Med(contamedia).Anisotropic[0].Sigma = this->ANIMATS.surfs(i).Sigma;
         sgg.Med(contamedia).Anisotropic[0].Mur = this->ANIMATS.surfs(i).mu / Mu0;
         sgg.Med(contamedia).Anisotropic[0].SigmaM = this->ANIMATS.surfs(i).SigmaM;
         tama2 = (this->ANIMATS.surfs(i).n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->ANIMATS.surfs(i).c2P(j).XI;
            punto.XE = this->ANIMATS.surfs(i).c2P(j).XE;
            punto.YI = this->ANIMATS.surfs(i).c2P(j).YI;
            punto.YE = this->ANIMATS.surfs(i).c2P(j).YE;
            punto.ZI = this->ANIMATS.surfs(i).c2P(j).ZI;
            punto.ZE = this->ANIMATS.surfs(i).c2P(j).ZE;
            orientacion = this->ANIMATS.surfs(i).c2P(j).or;
            numertag = searchtag(tagtype, this->ANIMATS.surfs(i).c2P(j).tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
         }
         tama3 = (this->ANIMATS.surfs(i).n_c1P);
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->ANIMATS.surfs(i).c1P(j).XI;
            punto.XE = this->ANIMATS.surfs(i).c1P(j).XI;
            punto.YI = this->ANIMATS.surfs(i).c1P(j).YI;
            punto.YE = this->ANIMATS.surfs(i).c1P(j).YI;
            punto.ZI = this->ANIMATS.surfs(i).c1P(j).ZI;
            punto.ZE = this->ANIMATS.surfs(i).c1P(j).ZI;
            orientacion = this->ANIMATS.surfs(i).c1P(j).or;
            numertag = searchtag(tagtype, this->ANIMATS.surfs(i).c1P(j).tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_

iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, &
            & contamedia);
         }
      }
      // LINs
      tama = this->ANIMATS.nLINS;
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         allocate(sgg.Med[contamedia].Anisotropic[1]);
         sgg.Med[contamedia].Is.Anisotropic = true;
         sgg.Med[contamedia].Priority = prior_IL;
         sgg.Med[contamedia].Anisotropic[1].Epr = this->ANIMATS.lins[i].eps / Eps0;
         sgg.Med[contamedia].Anisotropic[1].Sigma = this->ANIMATS.lins[i].Sigma;
         sgg.Med[contamedia].Anisotropic[1].Mur = this->ANIMATS.lins[i].mu / Mu0;
         sgg.Med[contamedia].Anisotropic[1].SigmaM = this->ANIMATS.lins[i].SigmaM;
         tama2 = this->ANIMATS.lins[i].n_c2P;
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->ANIMATS.lins[i].c2P[j].XI;
            punto.XE = this->ANIMATS.lins[i].c2P[j].XE;
            punto.YI = this->ANIMATS.lins[i].c2P[j].YI;
            punto.YE = this->ANIMATS.lins[i].c2P[j].YE;
            punto.ZI = this->ANIMATS.lins[i].c2P[j].ZI;
            punto.ZE = this->ANIMATS.lins[i].c2P[j].ZE;
            orientacion = this->ANIMATS.lins[i].c2P[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->ANIMATS.lins[i].c2P[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
         }
         tama3 = this->ANIMATS.lins[i].n_c1P;
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->ANIMATS.lins[i].c1P[j].XI;
            punto.XE = this->ANIMATS.lins[i].c1P[j].XI;
            punto.YI = this->ANIMATS.lins[i].c1P[j].YI;
            punto.YE = this->ANIMATS.lins[i].c1P[j].YI;
            punto.ZI = this->ANIMATS.lins[i].c1P[j].ZI;
            punto.ZE = this->ANIMATS.lins[i].c1P[j].ZI;
            orientacion = this->ANIMATS.lins[i].c1P[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->ANIMATS.lins[i].c1P[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
         }
      }
      // frequency dependent materials
      // bodies
      //
      tama = this->FRQDEPMATS.nvols;
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         fdgeom = nullptr;
         fdgeom = &this->FRQDEPMATS.vols[i];
         asignadisper(fdgeom);
         // geometry
         // !!!!!!!!!

         tama2 = this->FRQDEPMATS.vols[i].n_C;
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->FRQDEPMATS.vols[i].C[j].XI;
            punto.XE = this->FRQDEPMATS.vols[i].C[j].XE;
            punto.YI = this->FRQDEPMATS.vols[i].C[j].YI;
            punto.YE = this->FRQDEPMATS.vols[i].C[j].YE;
            punto.ZI = this->FRQDEPMATS.vols[i].C[j].ZI;
            punto.ZE = this->FRQDEPMATS.vols[i].C[j].ZE;
            numertag = searchtag(tagtype, this->FRQDEPMATS.vols[i].C[j].tag);
            CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                           media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                           Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                           Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                           Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                           Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                           Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
         }
      }


      // SURFs
      tama = this->FRQDEPMATS.nsurfs;
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         fdgeom = nullptr;
         fdgeom = &this->FRQDEPMATS.surfs[i];
         asignadisper(fdgeom);

         tama2 = this->FRQDEPMATS.surfs[i].n_C;
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->FRQDEPMATS.surfs[i].C[j].XI;
            punto.XE = this->FRQDEPMATS.surfs[i].C[j].XE;
            punto.YI = this->FRQDEPMATS.surfs[i].C[j].YI;
            punto.YE = this->FRQDEPMATS.surfs[i].C[j].YE;
            punto.ZI = this->FRQDEPMATS.surfs[i].C[j].ZI;
            punto.ZE = this->FRQDEPMATS.surfs[i].C[j].ZE;
            orientacion = this->FRQDEPMATS.surfs[i].C[j].or;
            numertag = searchtag(tagtype, this->FRQDEPMATS.surfs[i].C[j].tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
         }
      }
      // LINs
      tama = this->FRQDEPMATS.nLINS;
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         fdgeom = nullptr;
         fdgeom = &this->FRQDEPMATS.lins[i];
         asignadisper(fdgeom);

         tama2 = this->FRQDEPMATS.lins[i].n_C;
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->FRQDEPMATS.lins[i].C[j].XI;
            punto.XE = this->FRQDEPMATS.lins[i].C[j].XE;
            punto.YI = this->FRQDEPMATS.lins[i].C[j].YI;
            punto.YE = this->FRQDEPMATS.lins[i].C[j].YE;
            punto.ZI = this->FRQDEPMATS.lins[i].C[j].ZI;
            punto.ZE = this->FRQDEPMATS.lins[i].C[j].ZE;
            orientacion = this->FRQDEPMATS.lins[i].C[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->FRQDEPMATS.lins[i].C[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
         }
      }
      //
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // ISOTROPIC Multiports
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      inicontamedia = contamedia + 1;
      maxcontamedia = contamedia;
      tama = this->LossyThinSurfs.length;
      for (j = 1; j <= tama; ++j) {
         // carbon Multiports a guevo
         if (this->LossyThinSurfs.cs[j].numcapas == 0) {
            this->LossyThinSurfs.cs[j].numcapas = 1;
            allocate(this->LossyThinSurfs.cs[j].SigmaM[1]);
            allocate(this->LossyThinSurfs.cs[j].Sigma[1]);
            allocate(this->LossyThinSurfs.cs[j].EPS[1]);
            allocate(this->LossyThinSurfs.cs[j].MU[1]);
            allocate(this->LossyThinSurfs.cs[j].thk[1]);

            // _for_devia 090519
            allocate(this->LossyThinSurfs.cs[j].SigmaM_devia[1]);
            allocate(this->LossyThinSurfs.cs[j].Sigma_devia[1]);
            allocate(this->LossyThinSurfs.cs[j].EPS_devia[1]);
            allocate(this->LossyThinSurfs.cs[j].MU_devia[1]);
            allocate(this->LossyThinSurfs.cs[j].thk_devia[1]);
            // !!!

            this->LossyThinSurfs.cs[J].SigmaM[1] = 0.0_RKIND; // TRUCO PARA QUE CUANDO NO TENGA CAPAS (LECTURA DESDE FICHERO DE POLOS /RESIDUOS) NO PETE
            this->LossyThinSurfs.cs[J].Sigma[1] = 0.0_RKIND;
            this->LossyThinSurfs.cs[J].EPS[1] = EPS0;
            this->LossyThinSurfs.cs[J].MU[1] = MU0;
            this->LossyThinSurfs.cs[J].thk[1] = -1.0;

            // _for_devia 090519
            this->LossyThinSurfs.cs[J].SigmaM_devia[1] = 0.0_RKIND;
            this->LossyThinSurfs.cs[J].Sigma_devia[1] = 0.0_RKIND;
            this->LossyThinSurfs.cs[J].EPS_devia[1] = 0.0_RKIND;
            this->LossyThinSurfs.cs[J].MU_devia[1] = 0.0_RKIND;
            this->LossyThinSurfs.cs[J].thk_devia[1] = 0.0_RKIND;
            // !!!
            // !!!comentado el 120219 pq no se lleva bien con Semba !no entiendo ahora el comentario de malonyedispersive!!!120219
            // !!!     write(buff, '(a)')    'pre1_Error:  SGBC materials must have at least one layyer even in dummy for malonyedispersive'
            // !!!     call WarnErrReport (buff,.true.)
         }


         if (abs(this->LossyThinSurfs.cs[j].SigmaM[1]) <= 1.0e-2_RKIND) { // !!!ojoooo a 210319 manda guevos que tengamos que estar con el flag de la conductidad magnetica para llamar a SGBC todavia en 2015!!!
            this->LossyThinSurfs.cs[j].SigmaM = 0.0_RKIND;
            if (!mibc) {
               // if (this->LossyThinSurfs.cs(j)%numcapas >1) then
               //    write(buff, '(a)')    'pre1_Warning:  SGBC materials are just averaged for multilayered structures.'// &
               //    ' Use preferably -mibc instead.'
               //    call WarnErrReport (buff)
               // end if
               SGBC = true; // si la conductividad es 0.0_RKIND (o casi) utiliza directamente SGBC
               mibc = false;
            }
         }
         if (this->LossyThinSurfs.cs[j].SigmaM[1] >= 0.0_RKIND) {
            // SURFs (siempre son surfs)
            //
            tama2 = this->LossyThinSurfs.cs[j].nc;
            mincontamedia = maxcontamedia + 1;
            MultiportFile = trim(adjustl(this->LossyThinSurfs.cs[j].files)) + "_z11.txt";
            for (i = 1; i <= tama2; ++i) {
               orientacion = this->LossyThinSurfs.cs[j].C[i].or;
               punto.XI = this->LossyThinSurfs.cs[j].C[i].XI;
               punto.XE = this->LossyThinSurfs.cs[j].C[i].XE;
               punto.YI = this->LossyThinSurfs.cs[j].C[i].YI;
               punto.YE = this->LossyThinSurfs.cs[j].C[i].YE;
               punto.ZI = this->LossyThinSurfs.cs[j].C[i].ZI;
               punto.ZE = this->LossyThinSurfs.cs[j].C[i].ZE;
               existia = false;
               doexis: for (k = inicontamedia; k <= maxcontamedia; ++k) {
                  if (trim(adjustl(sgg.Med[k].multiport[1].multiportFileZ11)) == trim(adjustl(MultiportFile))) {
                     if (sgg.Med[k].multiport[1].Multiportdir == orientacion) {
                        contamedia = k;
                        existia = true;
                        break; // EXIT doexis
                     }
                  }
               }
               if (!existia) {
                  maxcontamedia = maxcontamedia + 1;
                  contamedia = maxcontamedia;
                  allocate(sgg.Med[contamedia].multiport[1]);
                  //
                  if ((this->LossyThinSurfs.cs[j].numcapas > 1) && SGBCDispersive) {
                     write(buff, *) "ERROR in SGBCs Number of layers >1 still unsupported for SGBCDispersive. ";
                     StopOnError(0, 0, buff);
                  }
                  //
                  allocate(sgg.Med[contamedia].Multiport[1].epr[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].mur[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].sigma[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].sigmam[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].width[1:this->LossyThinSurfs.cs[j].numcapas]);
                  // _for_devia 090519
                  allocate(sgg.Med[contamedia].Multiport[1].epr_devia[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].mur_devia[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].sigma_devia[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].sigmaM_devia[1:this->LossyThinSurfs.cs[j].numcapas],
                           sgg.Med[contamedia].Multiport[1].width_devia[1:this->LossyThinSurfs.cs[j].numcapas]);
                  // !!!
                  puntoXI = Max(punto.XI, Min(BoundingBox.XI, BoundingBox.XE));
                  puntoYI = Max(punto.YI, Min(BoundingBox.YI, BoundingBox.YE));
                  puntoZI = Max(punto.ZI, Min(BoundingBox.ZI, BoundingBox.ZE));

                  // !!!!!!!!estaba antes  maaaal. bug 140815verano
                  if (!((puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE))) {
                     puntoXI = sgg.allocDxI;
                     write(buff, *) "ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)", puntoXI, puntoYI, puntoZI, sgg.allocDxI, sgg.allocDyI, sgg.allocDzI;
                     WarnErrReport(buff, true);
                  }
                  if (!((puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE))) {
                     puntoYI = sgg.allocDyI;
                     write(buff, *) "ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)", puntoXI, puntoYI, puntoZI, sgg.allocDxI, sgg.allocDyI, sgg.allocDzI;
                     WarnErrReport(buff, true);
                  }
                  if (!((puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE))) {
                     puntoZI = sgg.allocDzI;
                     write(buff, *) "ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)", puntoXI, puntoYI, puntoZI, sgg.allocDxI, sgg.allocDyI, sgg.allocDzI;
                     WarnErrReport(buff, true);
                  }
                  dentro = (puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE) &&
                           (puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE) &&
                           (puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE);
                  delta = -1.0_RKIND;
                  if (DENTRO) {
                     switch (abs(this->LossyThinSurfs.cs[j].C[i].or)) {
                     case iEx:
                        delta = (sgg.DX(puntoXI) + sgg.DX(puntoXI - 1)) / 2.0_RKIND;
                        break;
                     case iEy:
                        delta = (sgg.DY(puntoYI) + sgg.Dy(puntoYI - 1)) / 2.0_RKIND;
                        break;
                     case iEz:
                        delta = (sgg.DZ(puntoZI) + sgg.Dz(puntoZI - 1)) / 2.0_RKIND;
                        break;
                     default:
                        write(buff, '(a)') "Buggy error 1 in preprocess composites. .";
                        STOPONERROR(layoutnumber, num_procs, buff);
                     }
                  } else {
                     write(buff, '(a)') "Buggy error 2 in preprocess composites. .";
                     STOPONERROR(layoutnumber, num_procs, buff);
                  }
                  sgg.Med[contamedia].Multiport[1].numcapas = this->LossyThinSurfs.cs[j].numcapas;
                  // el especificado
                  sgg.Med[contamedia].Multiport[1].Multiportdir = this->LossyThinSurfs.cs[j].C[i].or;
                  for (I_ = 1; I_ <= sgg.Med[contamedia].Multiport[1].numcapas; ++I_) {
                     if (sgg.Med[contamedia].Multiport[1].Multiportdir > 0) {
                        j_ = I_;
                     } else {
                        j_ = sgg.Med[contamedia].Multiport[1].numcapas - I_ + 1; // dale la vuelta (medios no simetricos) !0121
                     }
                     sgg.Med[contamedia].Multiport[1].epr[j_] = this->LossyThinSurfs.cs[j].eps[i_] / Eps0;
                     sgg.Med[contamedia].Multiport[1].mur[j_] = this->LossyThinSurfs.cs[j].mu[i_] / mu0;
                     sgg.Med[contamedia].Multiport[1].sigma[j_] = this->LossyThinSurfs.cs[j].Sigma[i_];
                     sgg.Med[contamedia].Multiport[1].sigmam[j_] = abs(this->LossyThinSurfs.cs[j].Sigmam[i_]);
                     sgg.Med[contamedia].Multiport[1].width[j_] = this->LossyThinSurfs.cs[j].thk[i_];

                     // _for_devia 090519
                     sgg.Med[contamedia].Multiport[1].epr_devia[j_] = this->LossyThinSurfs.cs[j].eps_devia[i_] / Eps0;
                     sgg.Med[contamedia].Multiport[1].mur_devia[j_] = this->LossyThinSurfs.cs[j].MU_devia[i_] / mu0;
                     sgg.Med[contamedia].Multiport[1].sigma_devia[j_] = this->LossyThinSurfs.cs[j].Sigma_devia[i_];
                     sgg.Med[contamedia].Multiport[1].sigmaM_devia[j_] = abs(this->LossyThinSurfs.cs[j].SigmaM_devia[i_]);
                     sgg.Med[contamedia].Multiport[1].width_devia[j_] = this->LossyThinSurfs.cs[j].thk_devia[i_];
                  }

                  rdummy = maxval(abs(this->LossyThinSurfs.cs[j].MU_devia)) + maxval(abs(this->LossyThinSurfs.cs[j].SigmaM_devia));
                  if (rdummy > 1.0e-15_RKIND) {
                     write(buff, '(a)') "Non null deviations found in sigmam or mu in composites. Still unsupported.";
                     STOPONERROR(layoutnumber, num_procs, buff);
                  }

// old pre 17/07/15
            // sgg.Med(contamedia).Multiport(1).transversalSpaceDelta = delta;
            sgg.Med(contamedia).Priority = prior_CS;
            sgg.Med(contamedia).Epr = this->LossyThinSurfs.cs(j).eps(1) / Eps0;
            sgg.Med(contamedia).Sigma = this->LossyThinSurfs.cs(j).Sigma(1);
            sgg.Med(contamedia).Mur = this->LossyThinSurfs.cs(j).mu(1) / Mu0;
            sgg.Med(contamedia).SigmaM = abs(this->LossyThinSurfs.cs(j).SigmaM(1)); // may be negative
            if (mibc) {
                sgg.Med(contamedia).Is.multiport = true;
                sgg.Med(contamedia).Is.Lossy = true;
            } else if (SGBC) {
                sgg.Med(contamedia).Is.SGBC = true;
                sgg.Med(contamedia).Is.Lossy = true;
                if (SGBCDispersive) sgg.Med(contamedia).Is.SGBCDispersive = true;
            } else {
                sprintf(buff, "Some -mibc -sgbc switch should be used for Composites.");
                STOPONERROR(layoutnumber, num_procs, buff);
            }
            sgg.Med(contamedia).Is.Dielectric = false;

            sgg.Med(contamedia).multiport(1).multiportFileZ11 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z11.txt";
            sgg.Med(contamedia).multiport(1).multiportFileZ22 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z22.txt";
            sgg.Med(contamedia).multiport(1).multiportFileZ12 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z12.txt";
            sgg.Med(contamedia).multiport(1).multiportFileZ21 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z12.txt";

            if (mibc) { // se trataran con MIBC !!!151161
                sgg.Med(contamedia).Is.SGBC = false;
                sgg.Med(contamedia).Is.SGBCDispersive = false;
                sgg.Med(contamedia).Is.Lossy = true;

                errnofile1 = false;
                errnofile2 = false;
                errnofile3 = false;
                errnofile4 = false;
                pozi = index(sgg.Med(contamedia).Multiport(1).multiportFileZ11, "_z11.txt");
                multiportfile2 = trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ11.substr(0, pozi - 1)));
                inquire(file=trim(adjustl(multiportfile2)), EXIST=errnofile1);
                inquire(file=trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ11)), EXIST=errnofile2);
                inquire(file=trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ12)), EXIST=errnofile3);
                inquire(file=trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ22)), EXIST=errnofile4);

                if (!(errnofile1 || (errnofile2 && errnofile3 && errnofile4))) {
                    sprintf(buff, "Neither New nor Old style mibc FILE %s EXISTS.", trim(adjustl(multiportfile2)));
                    WarnErrReport(buff, true);
                }
            }

            // end 09/07/13

            numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                            Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                            Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                            Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, contamedia);
        }
    }
    // Multiports
}
contamedia = maxcontamedia;
// end ISOTROPIC multiports

// ANISOTROPIC Multiports
inicontamedia = contamedia + 1;
maxcontamedia = contamedia;
tama = this->LossyThinSurfs.length;
for (j = 1; j <= tama; j++) {
    // carbon Multiports a guevo
    if (this->LossyThinSurfs.cs(j).SigmaM(1) < 0.0_RKIND) {
        // SURFs (siempre son surfs)
        tama2 = this->LossyThinSurfs.cs(j).nc;
        mincontamedia = maxcontamedia + 1;
        MultiportFile = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z11.txt";
        for (i = 1; i <= tama2; i++) {
            orientacion = this->LossyThinSurfs.cs(j).C(i).or;
            punto.XI = this->LossyThinSurfs.cs(j).C(i).XI;
            punto.XE = this->LossyThinSurfs.cs(j).C(i).XE;
            punto.YI = this->LossyThinSurfs.cs(j).C(i).YI;
            punto.YE = this->LossyThinSurfs.cs(j).C(i).YE;
            punto.ZI = this->LossyThinSurfs.cs(j).C(i).ZI;
            punto.ZE = this->LossyThinSurfs.cs(j).C(i).ZE;
            existia = false;
            for (k = inicontamedia; k <= maxcontamedia; k++) {
                if (trim(adjustl(sgg.Med(k).AnisMultiport(1).multiportFileZ11)) == trim(adjustl(MultiportFile))) {
                    if (sgg.Med(k).AnisMultiport(1).Multiportdir == orientacion) {
                        contamedia = k;
                        existia = true;
                        break;
                    }
                }
            }
            if (!existia) {
                maxcontamedia++;
                contamedia = maxcontamedia;
                allocate(sgg.Med(contamedia).AnisMultiport(1));

                if (this->LossyThinSurfs.cs(j).numcapas > 1) {
                    sprintf(buff, "pre1_ERROR:  Anisotropic multiport materials unsupported for multilayered structures.");
                    WarnErrReport(buff, true);
                }
                puntoXI = Max(punto.XI, Min(BoundingBox.XI, BoundingBox.XE)); // copiado de healer
                puntoYI = Max(punto.YI, Min(BoundingBox.YI, BoundingBox.YE));
                puntoZI = Max(punto.ZI, Min(BoundingBox.ZI, BoundingBox.ZE));

                // estaba antes maaaal. bug 140815verano
                if (!((puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE))) puntoXI = sgg.allocDxI;
                if (!((puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE))) puntoYI = sgg.allocDyI;
                if (!((puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE))) puntoZI = sgg.allocDzI;
                sprintf(buff, "ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)");
                WarnErrReport(buff, true);

                dentro = (puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE) &&
                         (puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE) &&
                         (puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE);
                delta = -1.0_RKIND;
                if (DENTRO) {
                    switch (abs(this->LossyThinSurfs.cs(j).C(i).or)) {
                        case iEx:
                            delta = (sgg.DX(puntoXI) + sgg.DX(puntoXI - 1)) / 2.0_RKIND;
                            break;
                        case iEy:
                            delta = (sgg.DY(puntoYI) + sgg.Dy(puntoYI - 1)) / 2.0_RKIND;
                            break;
                        case iEz:
                            delta = (sgg.DZ(puntoZI) + sgg.Dz(puntoZI - 1)) / 2.0_RKIND;
                            break;
                        default:
                            sprintf(buff, "Buggy error 1 in preprocess composites. .");
                            STOPONERROR(layoutnumber, num_procs, buff);
                    }
                } else {
                    sprintf(buff, "Buggy error 2 in preprocess composites. .");
                    STOPONERROR(layoutnumber, num_procs, buff);
                }
                sgg.Med(contamedia).AnisMultiport(1).Multiportdir = this->LossyThinSurfs.cs(j).C(i).or;
                sgg.Med(contamedia).AnisMultiport(1).epr = this->LossyThinSurfs.cs(j).eps / Eps0;
                sgg.Med(contamedia).AnisMultiport(1).mur = this->LossyThinSurfs.cs(j).mu / mu0;
                sgg.Med(contamedia).AnisMultiport(1).sigma = this->LossyThinSurfs.cs(j).Sigma;
                sgg.Med(contamedia).AnisMultiport(1).sigmam = abs(this->LossyThinSurfs.cs(j).Sigmam);
                sgg.Med(contamedia).AnisMultiport(1).width = this->LossyThinSurfs.cs(j).thk;
                sgg.Med(contamedia).Priority = prior_CS;
                sgg.Med(contamedia).Epr = this->LossyThinSurfs.cs(j).eps(1) / Eps0;
                sgg.Med(contamedia).Sigma = this->LossyThinSurfs.cs(j).Sigma(1);
                sgg.Med(contamedia).Mur = this->LossyThinSurfs.cs(j).mu(1) / Mu0;
                sgg.Med(contamedia).SigmaM = abs(this->LossyThinSurfs.cs(j).SigmaM(1)); // may be negative

                if (mibc) {
                    sgg.Med(contamedia).Is.Anismultiport = true;
                    sgg.Med(contamedia).Is.Lossy = true;
                } else {
                    sprintf(buff, "Some -mibc -sgbc switch should be used for Anisotropic Composites.");
                    STOPONERROR(layoutnumber, num_procs, buff);
                }

                sgg.Med(contamedia).Is.Dielectric = false;

                sgg.Med(contamedia).AnisMultiport(1).multiportFileZ11 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z11.txt";
                sgg.Med(contamedia).AnisMultiport(1).multiportFileZ22 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z22.txt";
                sgg.Med(contamedia).AnisMultiport(1).multiportFileZ12 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z12.txt";
                sgg.Med(contamedia).AnisMultiport(1).multiportFileZ21 = trim(adjustl(this->LossyThinSurfs.cs(j).files)) + "_z12.txt";
            }

            numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                            Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                            Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                            Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, contamedia);
        }
    }
    // Multiports
}
contamedia = maxcontamedia;
// end ANISOTROPIC multiports

// Multiports lossy padding
if (abs(attfactor - 1.0_RKIND) > 1.0e-12_RKIND) {
    tama = this->LossyThinSurfs.length;
    for (j = 1; j <= tama; j++) {
        tama2 = this->LossyThinSurfs.cs(j).nc;
        contamedia++;
        for (i = 1; i <= tama2; i++) {
            orientacion = this->LossyThinSurfs.cs(j).C(i).or;
            // ES UN free-space multiportpadding CON LA PRIORIDAD DE UN MULTIPORT con conductividad magnetica que luego se desanulara
            sgg.Med(contamedia).Priority = prior_CS;
            sgg.Med(contamedia).Is.multiport = false;
            sgg.Med(contamedia).Is.ANISmultiport = false;
            sgg.Med(contamedia).Is.MultiportPadding = true;
            sgg.Med(contamedia).Is.Lossy = true;
            sgg.Med(contamedia).Is.Dielectric = false;
            sgg.Med(contamedia).Epr = 1.0_RKIND; // this->LossyThinSurfs.cs(j).eps / Eps0
            sgg.Med(contamedia).Sigma = 0.0_RKIND; // abs(this->LossyThinSurfs.cs(j).Sigma) !may be negative
            sgg.Med(contamedia).Mur = 1.0_RKIND; // this->LossyThinSurfs.cs(j).mu / Mu0
            sgg.Med(contamedia).Is.Dielectric = true;
            // provisionalmente (luego se retocara con el sigmam correcto)
            sgg.Med(contamedia).SigmaM = 0.0_RKIND; // abs(this->LossyThinSurfs.cs(j).Sigma) !may be negative

            punto.XI = this->LossyThinSurfs.cs(j).C(i).XI;
            punto.XE = this->LossyThinSurfs.cs(j).C(i).XE;
            punto.YI = this->LossyThinSurfs.cs(j).C(i).YI;
            punto.YE = this->LossyThinSurfs.cs(j).C(i).YE;
            punto.ZI = this->LossyThinSurfs.cs(j).C(i).ZI;
            punto.ZE = this->LossyThinSurfs.cs(j).C(i).ZE;

            switch (Abs(orientacion)) {
                case iEx:
                    punto.XI = this->LossyThinSurfs.cs(j).C(i).XI;
                    punto.XE = this->LossyThinSurfs.cs(j).C(i).XI;
                    numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                    CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                                          media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                                          Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                                          Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                                          Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                                          Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, contamedia);
                    punto.XI = this->LossyThinSurfs.cs(j).C(i).XI - 1;
                    punto.XE = this->LossyThinSurfs.cs(j).C(i).XI - 1;
                    numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                    CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                                          media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                                          Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                                          Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                                          Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                                          Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, contamedia);
                    break;
                case iEy:
                    punto.YI = this->LossyThinSurfs.cs(j).C(i).YI;
                    punto.YE = this->LossyThinSurfs.cs(j).C(i).YI;
                    numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                    CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                                          media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                                          Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                                          Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                                          Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                                          Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, contamedia);
                    punto.YI = this->LossyThinSurfs.cs(j).C(i).YI - 1;
                    punto.YE = this->LossyThinSurfs.cs(j).C(i).YI - 1;
                    numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                    CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                                          media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                                          Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                                          Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                                          Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                                          Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, contamedia);
                    break;
                case iEz:
                    punto.ZI = this->LossyThinSurfs.cs(j).C(i).ZI;
                    punto.ZE = this->LossyThinSurfs.cs(j).C(i).ZI;
                    numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                    CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                                          media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE);

, Alloc_iEy_XI,
                     Alloc_iEy_XE, Alloc_iEy_YI,
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                     Alloc_iEz_YE, Alloc_iEz_ZI,
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,
                     Alloc_iHx_ZE, Alloc_iHy_XI,
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,
                     Alloc_iHz_XE, Alloc_iHz_YI,
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared,
                     BoundingBox, punto, orientacion,
                  & contamedia);
                  punto.ZI = this->LossyThinSurfs.cs(j).C(i).ZI - 1;
                  punto.ZE = this->LossyThinSurfs.cs(j).C(i).ZI - 1;
                  numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                  call CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy,
                     media.sggMiEz,
                  & media.sggMiHx, media.sggMiHy, media.sggMiHz,
                     Alloc_iEx_XI,
                  & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,
                     Alloc_iEy_XE, Alloc_iEy_YI,
                  & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                     Alloc_iEz_YE, Alloc_iEz_ZI,
                  & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,
                     Alloc_iHx_ZE, Alloc_iHy_XI,
                  & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,
                     Alloc_iHz_XE, Alloc_iHz_YI,
                  & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared,
                     BoundingBox, punto, orientacion,
                  & contamedia);
               } end select;
            } end do;
         } end do;
      } end if;
      //
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //wires
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      oldcontamedia = contamedIA;
      tama = this->twires.n_tw;
      do j = 1, tama
         contamedia = contamedia + 1;
        allocate(sgg.Med(contamedia).wire(1));
         sgg.Med(contamedia).Priority = prior_TW;

         //
         //background
         //
         sgg.Med(contamedia).Epr = sgg.Med(1).Epr;
         sgg.Med(contamedia).Sigma = sgg.Med(1).Sigma;
         sgg.Med(contamedia).Mur = sgg.Med(1).Mur;
         sgg.Med(contamedia).SigmaM = sgg.Med(1).SigmaM;
         sgg.Med(contamedia).Is.ThinWire = .TRUE.;
         sgg.Med(contamedia).Is.Dielectric = .FALSE.;
         sgg.Med(contamedia).wire(1).radius = this->twires.TW(j).RAD;
         sgg.Med(contamedia).wire(1).radius_devia = this->twires.TW(j).RAD_devia;
         if (boundwireradius) then
            if (sgg.Med(contamedia).wire(1).radius > maxwireradius) sgg.Med(contamedia).wire(1).radius = maxwireradius;
         end if;
         sgg.Med(contamedia).wire(1).R = this->twires.TW(j).RES;
         sgg.Med(contamedia).wire(1).l = this->twires.TW(j).IND;
         sgg.Med(contamedia).wire(1).C = this->twires.TW(j).CAP;
         sgg.Med(contamedia).wire(1).P_R = this->twires.TW(j).P_RES;
         sgg.Med(contamedia).wire(1).P_l = this->twires.TW(j).P_IND;

         sgg.Med(contamedia).wire(1).P_C = this->twires.TW(j).P_CAP;
         if (this->twires.TW(j).disp) then
            allocate(sgg.Med(contamedia).wire(1).disp(1));
            call asignawiredisper(sgg.Med(contamedia).wire(1).disp(1),
               this->twires.TW(j).dispfile);
         end if;
         sgg.Med(contamedia).wire(1).LeftEnd = this->twires.TW(j).LeftEnd;
         sgg.Med(contamedia).wire(1).RightEnd = this->twires.TW(j).RightEnd;
         sgg.Med(contamedia).wire(1).VsourceExists = .FALSE.;
         sgg.Med(contamedia).wire(1).IsourceExists = .FALSE.;
         //
         sgg.Med(contamedia).wire(1).HasAbsorbing_RightEnd = .FALSE.;
         sgg.Med(contamedia).wire(1).HasAbsorbing_LeftEnd = .FALSE.;
         sgg.Med(contamedia).wire(1).HasParallel_RightEnd = .FALSE.;
         sgg.Med(contamedia).wire(1).HasParallel_LeftEnd = .FALSE.;
         sgg.Med(contamedia).wire(1).HasSeries_RightEnd = .FALSE.;
         sgg.Med(contamedia).wire(1).HasSeries_LeftEnd = .FALSE.;
         sgg.Med(contamedia).wire(1).Parallel_R_RightEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_R_LeftEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_R_RightEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_R_LeftEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_L_RightEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_L_LeftEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_L_RightEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_L_LeftEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_C_RightEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_C_LeftEnd = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_C_RightEnd = 2.0e7_RKIND; //en corto 14/2/14
         sgg.Med(contamedia).wire(1).Series_C_LeftEnd = 2.0e7_RKIND; //en corto 14/2/14
//stoch
         sgg.Med(contamedia).wire(1).R_devia = this->twires.TW(j).RES_devia;
         sgg.Med(contamedia).wire(1).l_devia = this->twires.TW(j).IND_devia;
         sgg.Med(contamedia).wire(1).C_devia = this->twires.TW(j).CAP_devia;

         sgg.Med(contamedia).wire(1).Parallel_R_RightEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_R_LeftEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_L_RightEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_L_LeftEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_C_RightEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Parallel_C_LeftEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_R_RightEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_R_LeftEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_L_RightEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_L_LeftEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_C_RightEnd_devia = 0.0_RKIND;
         sgg.Med(contamedia).wire(1).Series_C_LeftEnd_devia = 0.0_RKIND;
//fin stoch
         //
         if (this->twires.TW(j).TL == MATERIAL_absorbing) then
            sgg.Med(contamedia).wire(1).HasAbsorbing_LeftEnd = .TRUE.;
         elseIF (this->twires.TW(j).TL == Parallel_CONS) then
            sgg.Med(contamedia).wire(1).HasParallel_LeftEnd = .TRUE.;
            sgg.Med(contamedia).wire(1).Parallel_R_LeftEnd = this->twires.TW(j).R_LeftEnd;
            sgg.Med(contamedia).wire(1).Parallel_L_LeftEnd = this->twires.TW(j).L_LeftEnd;
            sgg.Med(contamedia).wire(1).Parallel_C_LeftEnd = this->twires.TW(j).C_LeftEnd;
//
            sgg.Med(contamedia).wire(1).Parallel_R_LeftEnd_devia = this->twires.TW(j).R_LeftEnd_devia;
            sgg.Med(contamedia).wire(1).Parallel_L_LeftEnd_devia = this->twires.TW(j).L_LeftEnd_devia;
            sgg.Med(contamedia).wire(1).Parallel_C_LeftEnd_devia = this->twires.TW(j).C_LeftEnd_devia;

         ELSE if (this->twires.TW(j).TL == SERIES_CONS) then
            sgg.Med(contamedia).wire(1).HasSeries_LeftEnd = .TRUE.;
            sgg.Med(contamedia).wire(1).Series_R_LeftEnd = this->twires.TW(j).R_LeftEnd;
            sgg.Med(contamedia).wire(1).Series_L_LeftEnd = this->twires.TW(j).L_LeftEnd;
            sgg.Med(contamedia).wire(1).Series_C_LeftEnd = this->twires.TW(j).C_LeftEnd;
//
            sgg.Med(contamedia).wire(1).Series_R_LeftEnd_devia = this->twires.TW(j).R_LeftEnd_devia;
            sgg.Med(contamedia).wire(1).Series_L_LeftEnd_devia = this->twires.TW(j).L_LeftEnd_devia;
            sgg.Med(contamedia).wire(1).Series_C_LeftEnd_devia = this->twires.TW(j).C_LeftEnd_devia;
         ELSE if (this->twires.TW(j).TL == DISPERSIVE_CONS) then
            allocate(sgg.Med(contamedia).wire(1).disp_LeftEnd(1));
            call asignawiredisper(sgg.Med(contamedia).wire(1).disp_LeftEnd(1),
               this->twires.TW(j).dispfile_LeftEnd);
         end if;
         //

         if (this->twires.TW(j).TR == MATERIAL_absorbing) then
            sgg.Med(contamedia).wire(1).HasAbsorbing_RightEnd = .TRUE.;
         elseIF (this->twires.TW(j).TR == Parallel_CONS) then
            sgg.Med(contamedia).wire(1).HasParallel_RightEnd = .TRUE.;
            sgg.Med(contamedia).wire(1).Parallel_R_RightEnd = this->twires.TW(j).R_RightEnd;
            sgg.Med(contamedia).wire(1).Parallel_L_RightEnd = this->twires.TW(j).L_RightEnd;
            sgg.Med(contamedia).wire(1).Parallel_C_RightEnd = this->twires.TW(j).C_RightEnd;
//
            sgg.Med(contamedia).wire(1).Parallel_R_RightEnd_devia = this->twires.TW(j).R_RightEnd_devia;
            sgg.Med(contamedia).wire(1).Parallel_L_RightEnd_devia = this->twires.TW(j).L_RightEnd_devia;
            sgg.Med(contamedia).wire(1).Parallel_C_RightEnd_devia = this->twires.TW(j).C_RightEnd_devia;
         ELSE if (this->twires.TW(j).TR == SERIES_CONS) then
            sgg.Med(contamedia).wire(1).HasSeries_RightEnd = .TRUE.;
            sgg.Med(contamedia).wire(1).Series_R_RightEnd = this->twires.TW(j).R_RightEnd;
            sgg.Med(contamedia).wire(1).Series_L_RightEnd = this->twires.TW(j).L_RightEnd;
            sgg.Med(contamedia).wire(1).Series_C_RightEnd = this->twires.TW(j).C_RightEnd;
//
            sgg.Med(contamedia).wire(1).Series_R_RightEnd_devia = this->twires.TW(j).R_RightEnd_devia;
            sgg.Med(contamedia).wire(1).Series_L_RightEnd_devia = this->twires.TW(j).L_RightEnd_devia;
            sgg.Med(contamedia).wire(1).Series_C_RightEnd_devia = this->twires.TW(j).C_RightEnd_devia;

         ELSE if (this->twires.TW(j).TR == DISPERSIVE_CONS) then
            allocate(sgg.Med(contamedia).wire(1).disp_RightEnd(1));
            call asignawiredisper(sgg.Med(contamedia).wire(1).disp_RightEnd(1),
               this->twires.TW(j).dispfile_RightEnd);
         end if;

//stoch
         rdummy = abs(sgg.Med(contamedia).wire(1).radius_devia)
            + abs(sgg.Med(contamedia).wire(1).l_devia)
            + abs(sgg.Med(contamedia).wire(1).C_devia)
            + abs(sgg.Med(contamedia).wire(1).Parallel_L_RightEnd_devia)
            + abs(sgg.Med(contamedia).wire(1).Parallel_L_LeftEnd_devia)
            + abs(sgg.Med(contamedia).wire(1).Parallel_C_RightEnd_devia)
            + abs(sgg.Med(contamedia).wire(1).Parallel_C_LeftEnd_devia)
            + abs(sgg.Med(contamedia).wire(1).Series_L_RightEnd_devia)
            + abs(sgg.Med(contamedia).wire(1).Series_L_LeftEnd_devia)
            + abs(sgg.Med(contamedia).wire(1).Series_C_RightEnd_devia)
            + abs(sgg.Med(contamedia).wire(1).Series_C_LeftEnd_devia);
         if (rdummy > 1.0e-15_RKIND) then
            write(buff, '(a)') 'Non null deviations found in L, C or radius in wires stoch. Still unsupported.';
            call STOPONERROR(layoutnumber, num_procs, buff);
         end if;
//fin stoch
         //
         //esto se soportaba desde versiones antiguas (hilos de un solo segmento. Por error se descomento en la R2417 cuando se trabajo en lo del strictnfde tras vuelta de madrid
         //vuelvo a comentarlo porque si que tenemos la capacidad de hilos de un solo segmento
         //
         //        tama2 = this->twires.TW(j).N_TWC;
         //        if (tama2 == 1) then
         //           call stoponerror(layoutnumber, num_procs, 'A WIRE must have at least two segments');
         //        end if;
         //
         //esto no es ya necesario porque lo calculo yo luego en el wires
         //!record the LeftEnd and RightEnd coordinates (first and last points)
         //!
         //        sgg.Med(contamedia).wire(1).LextremoI = this->twires.TW(j).TWC(1).i;
         //        sgg.Med(contamedia).wire(1).LextremoJ = this->twires.TW(j).TWC(1).j;
         //        sgg.Med(contamedia).wire(1).LextremoK = this->twires.TW(j).TWC(1).k;
         //        orientacionL = this->twires.TW(j).TWC(1).D;
         //        sgg.Med(contamedia).wire(1).RextremoI = this->twires.TW(j).TWC(tama2).i;
         //        sgg.Med(contamedia).wire(1).RextremoJ = this->twires.TW(j).TWC(tama2).j;
         //        sgg.Med(contamedia).wire(1).RextremoK = this->twires.TW(j).TWC(tama2).k;
         //        orientacionR = this->twires.TW(j).TWC(tama2).D;
         //        //
         //!correct each ending
         //        numminus = 0;
         //        do i = 2, tama2 - 1 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         //            punto.XI = this->twires.TW(j).TWC(i).i;
         //            punto.YI = this->twires.TW(j).TWC(i).j;
         //            punto.ZI = this->twires.TW(j).TWC(i).k;
         //            orientacion = this->twires.TW(j).TWC(i).D;
         //            select case (orientacion)
         //            case (iEx)
         //                if ((punto.YI == sgg.Med(contamedia).wire(1).LextremoJ) .and.
         //                          (punto.ZI == sgg.Med(contamedia).wire(1).LextremoK)) then
         //                    if ((orientacion /= orientacionL) .and. (punto.XI == sgg.Med(contamedia).wire(1).LextremoI)) numminus = numminus + 1; !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         //                    if (punto.XI + 1 == sgg.Med(contamedia).wire(1).LextremoI) numminus = numminus + 1;
         //                end if;
         //            case (iEy)
         //                if ((punto.XI == sgg.Med(contamedia).wire(1).LextremoI) .and.
         //                          (punto.ZI == sgg.Med(contamedia).wire(1).LextremoK)) then
         //                    if ((orientacion /= orientacionL) .and. (punto.YI == sgg.Med(contamedia).wire(1).LextremoJ)) numminus = numminus + 1;
         //                    if (punto.YI + 1 == sgg.Med(contamedia).wire(1).LextremoJ) numminus = numminus + 1;
         //                end if;
         //            case (iEz)
         //                if ((punto.YI == sgg.Med(contamedia).wire(1).LextremoJ) .and.
         //                          (punto.XI == sgg.Med(contamedia).wire(1).LextremoI)) then
         //                    if ((orientacion /= orientacionL) .and. (punto.ZI == sgg.Med(contamedia).wire(1).LextremoK)) numminus = numminus + 1;
         //                    if (punto.ZI + 1 == sgg.Med(contamedia).wire(1).LextremoK) numminus = numminus + 1;
         //                end if;
         //            end select;
         //
         //        end do;
         //        if (numminus >= 1) then !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         //              select case (this->twires.TW(j).TWC(1).D)
         //              case (iEx) !si son iguales a 2 es cerrado
         //                  sgg.Med(contamedia).wire(1).LextremoI = sgg.Med(contamedia).wire(1).LextremoI + 1;
         //              case (iEy)
         //                  sgg.Med(contamedia).wire(1).LextremoJ = sgg.Med(contamedia).wire(1).LextremoJ + 1;
         //              case (iEz)
         //                  sgg.Med(contamedia).wire(1).LextremoK = sgg.Med(contamedia).wire(1).LextremoK + 1;
         //              end select;
         //        end if;
         //        //
         //!correct each ending
         //        numminus = 0;
         //        do i = 2, tama2 - 1 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         //            punto.XI = this->twires.TW(j).TWC(i).i;
         //            punto.YI = this->twires.TW(j).TWC(i).j;
         //            punto.ZI = this->twires.TW(j).TWC(i).k;
         //            orientacion = this->twires.TW(j).TWC(i).D;
         //            select case (orientacion)
         //            case (iEx)
         //                if ((punto.YI == sgg.Med(contamedia).wire(1).RextremoJ) .and.
         //                          (punto.ZI == sgg.Med(contamedia).wire(1).RextremoK)) then
         //                    if ((orientacion /= orientacionR) .and. (punto.XI == sgg.Med(contamedia).wire(1).RextremoI)) numminus = numminus + 1;
         //                    if (punto.XI + 1 == sgg.Med(contamedia).wire(1).RextremoI) numminus = numminus + 1;
         //                end if;
         //            case (iEy)
         //                if ((punto.XI == sgg.Med(contamedia).wire(1).RextremoI) .and.
         //                          (punto.ZI == sgg.Med(contamedia).wire(1).RextremoK)) then
         //                    if ((orientacion /= orientacionR) .and. (punto.YI == sgg.Med(contamedia).wire(1).RextremoJ)) numminus = numminus + 1;
         //                    if (punto.YI + 1 == sgg.Med(contamedia).wire(1).RextremoJ) numminus = numminus + 1;
         //                end if;
         //            case (iEz)
         //                if ((punto.YI == sgg.Med(contamedia).wire(1).RextremoJ) .and.
         //                          (punto.XI == sgg.Med(contamedia).wire(1).RextremoI)) then
         //                    if ((orientacion /= orientacionR) .and. (punto.ZI == sgg.Med(contamedia).wire(1).RextremoK)) numminus = numminus + 1;
         //                    if (punto.ZI + 1 == sgg.Med(contamedia).wire(1).RextremoK) numminus = numminus + 1;
         //                end if;
         //            end select;
         //        end do;
         //        if ((numminus >= 1) .or. (tama2 == 1)) then  !bug ca295 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
         //              select case (this->twires.TW(j).TWC(tama2).D)
         //              case (iEx)
         //                  sgg.Med(contamedia).wire(1).RextremoI = sgg.Med(contamedia).wire(1).RextremoI + 1;
         //!si son iguales a 2 es cerrado
         //              case (iEy)
         //                  sgg.Med(contamedia).wire(1).RextremoJ = sgg.Med(contamedia).wire(1).RextremoJ + 1;
         //              case (iEz)
         //                  sgg.Med(contamedia).wire(1).RextremoK = sgg.Med(contamedia).wire(1).RextremoK + 1;
         //              end select;
         //        end if;
      } end do; !del tama

      //preanalisis de hilos embeddeds en materiales  antes de asignarlos
      tama = this->twires.n_tw;
      paraerrhilo = .false.;
      do j1 = 1, tama
         tama2 = this->twires.TW(j1).N_TWC;
         do i1 = 1, tama2
            i = this->twires.TW(j1).TWC(i1).i;
            j = this->twires.TW(j1).TWC(i1).j;
            k = this->twires.TW(j1).TWC(i1).k;
            orientacion = this->twires.TW(j1).TWC(i1).D;
            OrigIndex = this->twires.TW(j1).TWC(i1).nd;
            if ((i >= BoundingBox.XI) .AND. (i < BoundingBox.XE) .AND.
            &    (j >= BoundingBox.YI) .AND. (j < BoundingBox.YE) .AND.
            &    (k >= BoundingBox.ZI) .AND. (k < BoundingBox.ZE)) then
               if (i > BoundingBox.XI) then
                  imenos1 = i - 1;
               else
                  imenos1 = i;
               end if;
               if (j > BoundingBox.YI) then
                  jmenos1 = j - 1;
               else
                  jmenos1 = j;
               end if;
               if (k > BoundingBox.ZI) then
                  kmenos1 = k - 1;
               else
                  kmenos1 = k;
               end if;
               select c

case (iEx):
                  if ((media%sggMiEx(i,j,k) ==0).or.(sgg%med(media%sggMiEx(i,j,k) )%is%pec)) {
                     paraerrhilo=.true.;
                     write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   x-WIRE at ',OrigIndex, i, j, k,' embedded within PEC';
                     if (verbose) call WarnErrReport (buff);
                  } elseif (media%sggMiEx(i,j,k) /= 1) {
                     islossy = (sgg%Med(media%sggMiEx(i,j,k))%Sigma /= 0.0_RKIND);
                     if (islossy) {
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k, &
                           ' embedded within LOSSY medium ', &
                           media%sggMiEx(i,j,k);
                     } else {
                        write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k,' embedded within medium ', &
                           media%sggMiEx(i,j,k);
                     }
                     if (verbose) call WarnErrReport (buff);
                  }
                  if ((((media%sggMiEy(i  ,j,k) ==0).or.(sgg%med(media%sggMiEy(i  ,j,k) )%is%pec)).or. &
                     ((media%sggMiEz(i  ,j,k) ==0).or.(sgg%med(media%sggMiEz(i  ,j,k) )%is%pec)).or. &
                     ((media%sggMiEy(i  ,jmenos1,k) ==0).or.(sgg%med(media%sggMiEy(i  ,jmenos1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i  ,j,kmenos1) ==0).or.(sgg%med(media%sggMiEz(i  ,j,kmenos1) )%is%pec))).and. &
                  &     ((media%sggMiEx(i  ,j,k) /=0).and.(.not.(sgg%med(media%sggMiEx(i  ,j,k) )%is%pec)))) {
                     if ((i1 /= 1) .and. (i1 /= tama2)) { !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of x-WIRE at  ',OrigIndex, i, j, k, &
                           ' touching PEC';
                        if (verbose) call WarnErrReport (buff);
                     } else {
                        continue;
                        !write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC';
                        !if (verbose) call WarnErrReport (buff);
                     }
                  } elseif ((((media%sggMiEy(i+1,j,k) ==0).or.(sgg%med(media%sggMiEy(i+1,j,k) )%is%pec)).or.&
                     ((media%sggMiEz(i+1,j,k) ==0).or.(sgg%med(media%sggMiEz(i+1,j,k) )%is%pec)).or. &
                     ((media%sggMiEy(i+1,jmenos1,k) ==0).or.(sgg%med(media%sggMiEy(i+1,jmenos1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i+1,j,kmenos1) ==0).or.(sgg%med(media%sggMiEz(i+1,j,kmenos1) )%is%pec))).and. &
                  &         ((media%sggMiEx(i  ,j,k) /=0).and.(.not.(sgg%med(media%sggMiEx(i  ,j,k) )%is%pec)))) {
                     if ((i1 /= 1) .and. (i1 /= tama2)) { !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of x-WIRE at  ',OrigIndex, i+1, j, k, &
                           ' touching PEC';
                        if (verbose) call WarnErrReport (buff);
                     } else {
                        continue;
                        !write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i+1, j, k,' touching PEC';
                        !if (verbose) call WarnErrReport (buff);
                     }
                  } elseif (((media%sggMiEy(i  ,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i,j,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEz(i  ,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEz(i,j,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEy(i+1,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i+1, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i+1,j,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEz(i+1,j,k) /= 1)).and. &
                     (media%sggMiEx(i  ,j,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: x-WIRE at ',OrigIndex, i+1, j, k,' touching medium ', &
                     &                                  media%sggMiEz(i+1,j,k);
                     if (verbose) call WarnErrReport (buff);
                  }
                case (iEy):
                  if ((media%sggMiEy(i,j,k) ==0).or.(sgg%med(media%sggMiEy(i,j,k) )%is%pec)) {
                     paraerrhilo=.true.;
                     write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   y-WIRE at ',OrigIndex, i, j, k,' embedded within PEC';
                     if (verbose) call WarnErrReport (buff);
                  } elseif (media%sggMiEy(i,j,k) /= 1) {
                     islossy = (sgg%Med(media%sggMiEy(i,j,k))%Sigma /= 0.0_RKIND);
                     if (islossy) {
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: Y-WIRE at ',OrigIndex, i, j, k,' embedded within LOSSY medium ', &
                           media%sggMiEY(i,j,k);
                     } else {
                        write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j, k,' embedded within medium ', &
                           media%sggMiEy(i,j,k);
                     }
                     if (verbose) call WarnErrReport (buff);
                  }
                  if ((((media%sggMiEx(i,j  ,k) ==0).or.(sgg%med(media%sggMiEx(i,j  ,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j,k  ) ==0).or.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEx(imenos1,j  ,k) ==0).or.(sgg%med(media%sggMiEx(imenos1,j  ,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j,kmenos1  ) ==0).or.(sgg%med(media%sggMiEz(i,j,kmenos1  ) )%is%pec))).and. &
                  &     ((media%sggMiEy(i,j  ,k) /=0).and.(.not.(sgg%med(media%sggMiEy(i,j  ,k) )%is%pec)))) {
                     if ((i1 /= 1) .and. (i1 /= tama2)) { !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of y-WIRE at ',OrigIndex, i, j, k, &
                           ' touching PEC';
                        if (verbose) call WarnErrReport (buff);
                     } else {
                        continue;
                        !write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC';
                        !if (verbose) call WarnErrReport (buff);
                     }
                  } elseif ((((media%sggMiEx(i,j+1,k) ==0).or.(sgg%med(media%sggMiEx(i,j+1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j+1,k) ==0).or.(sgg%med(media%sggMiEz(i,j+1,k) )%is%pec)).or. &
                     ((media%sggMiEx(imenos1,j+1,k) ==0).or.(sgg%med(media%sggMiEx(imenos1,j+1,k) )%is%pec)).or. &
                     ((media%sggMiEz(i,j+1,kmenos1) ==0).or.(sgg%med(media%sggMiEz(i,j+1,kmenos1) )%is%pec))).and. &
                  &         ((media%sggMiEy(i,j  ,k) /=0).and.(.not.(sgg%med(media%sggMiEy(i,j  ,k) )%is%pec)))) {
                     if ((i1 /= 1) .and. (i1 /= tama2)) { !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of y-WIRE at ',OrigIndex, i, j+1, k, &
                           ' touching PEC';
                        if (verbose) call WarnErrReport (buff);
                     } else {
                        continue;
                        !write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j+1, k,' touching PEC';
                        !if (verbose) call WarnErrReport (buff);
                     }
                  } elseif (((media%sggMiEx(i,j  ,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEz(i,j  ,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEz(i,j,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEx(i,j+1,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j+1, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j+1,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEz(i,j+1,k) /= 1)).and. &
                  &         (media%sggMiEy(i,j  ,k) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: y-WIRE at ',OrigIndex, i, j+1, k,' touching medium ', &
                     &                                  media%sggMiEz(i,j+1,k);
                     if (verbose) call WarnErrReport (buff);
                  }
                case (iEz):
                  if ((media%sggMiEz(i,j,k) ==0).or.(sgg%med(media%sggMiEz(i,j,k) )%is%pec)) {
                     paraerrhilo=.true.;
                     write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   z-WIRE at ',OrigIndex, i, j, k,' embedded within PEC';
                     if (verbose) call WarnErrReport (buff);
                  } elseif (media%sggMiEz(i,j,k) /= 1) {
                     islossy = (sgg%Med(media%sggMiEz(i,j,k))%Sigma /= 0.0_RKIND);
                     if (islossy) {
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: Y-WIRE at ',OrigIndex, i, j, k, &
                           ' embedded within LOSSY medium ', &
                           media%sggMiEz(i,j,k);
                     } else {
                        write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' embedded within medium ', &
                           media%sggMiEz(i,j,k);
                     }
                     if (verbose) call WarnErrReport (buff);
                  }
                  if ((((media%sggMiEx(i,j,k  ) ==0).or.(sgg%med(media%sggMiEx(i,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEy(i,j,k  ) ==0).or.(sgg%med(media%sggMiEy(i,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEx(imenos1,j,k  ) ==0).or.(sgg%med(media%sggMiEx(imenos1,j,k  ) )%is%pec)).or. &
                     ((media%sggMiEy(i,jmenos1,k  ) ==0).or.(sgg%med(media%sggMiEy(i,jmenos1,k  ) )%is%pec))).and. &
                  &     ((media%sggMiEz(i,j,k  ) /=0).and.(.not.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec)))) {
                     if ((i1 /= 1) .and. (i1 /= tama2)) { !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of z-WIRE at ',OrigIndex, i, j, k, &
                           ' touching PEC';
                        if (verbose)  call WarnErrReport (buff);
                     } else {
                        continue;
                        !write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC';
                        !if (verbose) call WarnErrReport (buff);
                     }
                  } elseif ((((media%sggMiEx(i,j  ,k+1) ==0).or.(sgg%med(media%sggMiEx(i,j  ,k+1) )%is%pec)).or. &
                     ((media%sggMiEy(i,j  ,k+1) ==0).or.(sgg%med(media%sggMiEy(i,j  ,k+1) )%is%pec)).or.   &
                  &         ((media%sggMiEx(imenos1,j,k+1) ==0).or.(sgg%med(media%sggMiEx(imenos1,j,k+1) )%is%pec)).or. &
                     ((media%sggMiEy(i,jmenos1,k+1) ==0).or.(sgg%med(media%sggMiEy(i,jmenos1,k+1) )%is%pec))).and. &
                  &         (((media%sggMiEz(i,j,k  ) /=0).and.(.not.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec))).or.(.not.(sgg%med(media%sggMiEz(i,j,k  ) )%is%pec)))) {
                     if ((i1 /= 1) .and. (i1 /= tama2)) { !solo en LeftEnd y RightEnd pueden tocar
                        paraerrhilo=.true.;
                        write(buff, '(a,i7,3i5,a)')    'pre1_WARNING:   intermediate node of z-WIRE at ',OrigIndex, i, j, k+1, &
                           ' touching PEC';
                        if (verbose) call WarnErrReport (buff);
                     } else {
                        continue;
                        !write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k+1,' touching PEC';
                        !if (verbose) call WarnErrReport (buff);
                     }
                  } elseif (((media%sggMiEx(i,j,k  ) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEy(i,j,k  ) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i,j,k);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEy(i,j,k+1) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEy(i,j,k+1);
                     if (verbose) call WarnErrReport (buff);
                  } elseif (((media%sggMiEx(i,j,k+1) /= 1)).and. &
                  &         (media%sggMiEz(i,j,k  ) == 1)) {
                     write(buff, '(a,i7,3i5,a,i5)') 'pre1_WARNING: z-WIRE at ',OrigIndex, i, j, k,' touching medium ', &
                     &                                  media%sggMiEx(i,j,k+1);
                     if (verbose) call WarnErrReport (buff);
                  }
               end select
            end if
         end do
      end do
      !lo dejo que siga !luego el wires parara
      !if (paraerrhilo.and.(.not.groundwires)) then
      !    buff='Revise WIRE intersections!'
      !    call STOPONERROR(layoutnumber,num_procs,buff)
      !end if
      !end preanalisis
      !perform the assignments
      bboxwirxi=  2**20;
      bboxwirxe=-(2**20);
      bboxwirYi=  2**20;
      bboxwirYe=-(2**20);
      bboxwirZi=  2**20;
      bboxwirZe=-(2**20);
      contamedia = oldcontamedia;
      tama = this%twires%n_tw;
      for (j = 1; j <= tama; j++) {
         contamedia = contamedia + 1;
         tama2 = this%twires%TW(j)%N_TWC;
         TAMA2BIS=0;
         for (i = 1; i <= tama2; i++) {
            punto%XI = this%twires%TW(j)%TWC(i)%i;
            punto%XE = this%twires%TW(j)%TWC(i)%i;
            punto%YI = this%twires%TW(j)%TWC(i)%j;
            punto%YE = this%twires%TW(j)%TWC(i)%j;
            punto%ZI = this%twires%TW(j)%TWC(i)%k;
            punto%ZE = this%twires%TW(j)%TWC(i)%k;
            !!!!!!!!!los clipeo agresivamente si lo lanzo con -CLIPREGION para que no me den problema 06/07/15 (solo sirve para debugeo y con el -wiresflavor holland (old))
            if (((((punto%XI-2 >  SINPML_fullsize(iHx)%XI).and.(punto%XI+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YI-2 >  SINPML_fullsize(iHy)%YI).and.(punto%YI+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZI-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%zI+2 < SINPML_fullsize(iHz)%ZE).and. &
               (punto%XE-2 >  SINPML_fullsize(iHx)%XI).and.(punto%Xe+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YE-2 >  SINPML_fullsize(iHy)%YI).and.(punto%Ye+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZE-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%Ze+2 < SINPML_fullsize(iHz)%ZE))).or.(!CLIPREGION))) TAMA2BIS=TAMA2bis+1;
         }
         !
        allocate(sgg%Med(contamedia)%wire(1)%segm(1:TAMA2BIS));
         sgg%Med(contamedia)%wire(1)%numsegmentos = TAMA2BIS;
        allocate(sgg%Med(contamedia)%wire(1)%VSource(1:TAMA2BIS));
        allocate(sgg%Med(contamedia)%wire(1)%ISource(1:TAMA2BIS));
         CONTAVOLT=0;
         CONTACURR=0;

         TAMA2BIS=0;
         hilosbarre: do i = 1, tama2 {
            punto%XI = this%twires%TW(j)%TWC(i)%i;
            punto%XE = this%twires%TW(j)%TWC(i)%i;
            punto%YI = this%twires%TW(j)%TWC(i)%j;
            punto%YE = this%twires%TW(j)%TWC(i)%j;
            punto%ZI = this%twires%TW(j)%TWC(i)%k;
            punto%ZE = this%twires%TW(j)%TWC(i)%k;
!!!sgg250418
!!!bug 2018
!!bug que aparece cuando en hilos de dos segmentos con ambos identicos. Lo que hago es clipearlo directamente.
            if ((i==2).and.(tama2==2)) {
               if  (((this%twires%TW(j)%TWC(i)%i).eq.(this%twires%TW(j)%TWC(i-1)%i)).and. &
                  ((this%twires%TW(j)%TWC(i)%j).eq.(this%twires%TW(j)%TWC(i-1)%j)).and. &
                  ((this%twires%TW(j)%TWC(i)%k).eq.(this%twires%TW(j)%TWC(i-1)%k)).and. &
                  ((this%twires%TW(j)%TWC(i)%D).eq.(this%twires%TW(j)%TWC(i-1)%D))) {
                  sgg%Med(contamedia)%wire(1)%numsegmentos = 1;
                  write(buff, '(a,1i7,a)')    'pre1_SEVEREWARNING: removing a repeteated segment from a 2-segment wire. Index= ',this%twires%TW(j)%TWC(i)%nd,'. Double check that no wire probes was attached to it';
                  call WarnErrReport (buff);
                  exit hilosbarre;
               }
            }
!!!fin 250418

            !!!!!!!!!los clipeo agresivamente si lo lanzo con -CLIPREGION para que no me den problema 06/07/15 (solo sirve para debugeo y con el -wiresflavor holland (old))
            if (!( &
               (((punto%XI-2 >  SINPML_fullsize(iHx)%XI).and.(punto%XI+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YI-2 >  SINPML_fullsize(iHy)%YI).and.(punto%YI+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZI-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%zI+2 < SINPML_fullsize(iHz)%ZE).and. &
               (punto%XE-2 >  SINPML_fullsize(iHx)%XI).and.(punto%Xe+2 < SINPML_fullsize(iHx)%XE).and. &
               (punto%YE-2 >  SINPML_fullsize(iHy)%YI).and.(punto%Ye+2 < SINPML_fullsize(iHy)%YE).and. &
               (punto%ZE-2 >  SINPML_fullsize(iHz)%ZI).and.(punto%Ze+2 < SINPML_fullsize(iHz)%ZE))).or.(!CLIPREGION))) CYCLE hilosbarre;

            TAMA2BIS=TAMA2BIS+1;
            orientacion = this%twires%TW(j)%TWC(i)%D;
            origindex=    this%twires%TW(j)%TWC(i)%nd;
            !
            !
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%i = punto%XI;
            sgg%Med(contamedia)%wire(1)%SEGM(TAMA2BIS)%j = punto%YI;

sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->k = punto.ZI;

            // 2014 para informacion bbox hilos
            if (punto.XI < bboxwirxi) bboxwirXI = punto.XI;
            if (punto.XE > bboxwirxE) bboxwirXE = punto.XE;
            if (punto.YI < bboxwirYi) bboxwirYI = punto.YI;
            if (punto.YE > bboxwirYE) bboxwirYE = punto.YE;
            if (punto.ZI < bboxwirZi) bboxwirZI = punto.ZI;
            if (punto.ZE > bboxwirZE) bboxwirZE = punto.ZE;
            // !!!!!

            sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->ori = orientacion;
            sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->origindex = origindex;
            sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->Is_LeftEnd = false;
            sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->Is_RightEnd = false;
            sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->repetido = false; // luego el preprocesador del wires cambia esto
            if (i == 1) sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->Is_LeftEnd = true;
            if (i == tama2) sgg->Med[contamedia]->wire[0]->SEGM[TAMA2BIS]->Is_RightEnd = true;
            //
            isathinwire = true;
            numertag = searchtag(tagtype, this->twires->TW[j]->TWC[i]->tag);
            CreateLineMM(layoutnumber, media->sggMtag, tag_numbers, numertag, media->sggMiEx, media->sggMiEy, media->sggMiEz,
                media->sggMiHx, media->sggMiHy, media->sggMiHz, Alloc_iEx_XI,
                Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg->Med, sgg->NumMedia, sgg->EShared, BoundingBox, punto, orientacion,
                contamedia, isathinwire, verbose, numeroasignaciones);
            if ((trim(adjustl(this->twires->TW[j]->TWC[i]->SRCTYPE)) == F_SOURCE_VOLTAGE) ||
                (trim(adjustl(this->twires->TW[j]->TWC[i]->SRCTYPE)) == F_SOURCE_CURRENT)) {
                // !!!!!!!!!!!!!
                if (trim(adjustl(this->twires->TW[j]->TWC[i]->SRCTYPE)) == F_SOURCE_VOLTAGE) {
                    CONTAVOLT = CONTAVOLT + 1;
                    sgg->Med[contamedia]->wire[0]->VsourceExists = true;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->Multiplier = this->twires->TW[j]->TWC[i]->m;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->Resistance = 0.0_RKIND;
                    // not provided by .nfde but supported by the simulation though untested
                    //
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->fichero.name = trim(adjustl(this->twires->TW[j]->TWC[i]->SRCFILE));
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->i = this->twires->TW[j]->TWC[i]->i;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->j = this->twires->TW[j]->TWC[i]->j;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->k = this->twires->TW[j]->TWC[i]->k;
                } else if (trim(adjustl(this->twires->TW[j]->TWC[i]->SRCTYPE)) == F_SOURCE_CURRENT) {
                    CONTAVOLT = CONTAVOLT + 1;
                    sgg->Med[contamedia]->wire[0]->VsourceExists = true;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->Multiplier = 1.0e22_RKIND;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->Resistance = 1.0e22_RKIND;
                    // not provided by .nfde but supported by the simulation though untested
                    //
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->fichero.name = trim(adjustl(this->twires->TW[j]->TWC[i]->SRCFILE));
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->i = this->twires->TW[j]->TWC[i]->i;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->j = this->twires->TW[j]->TWC[i]->j;
                    sgg->Med[contamedia]->wire[0]->VSource[CONTAVOLT]->k = this->twires->TW[j]->TWC[i]->k;
                }
            } else if (trim(adjustl(this->twires->TW[j]->TWC[i]->SRCTYPE)) != "None") {
                sprintf(buff, "WRONG type of wire source %s", trim(adjustl(this->twires->TW[j]->TWC[i]->SRCTYPE)));
                stoponerror(layoutnumber, num_procs, buff);
            }
        } while (hilosbarre);
        sgg->Med[contamedia]->wire[0]->NUMVOLTAGESOURCES = CONTAVOLT;
        sgg->Med[contamedia]->wire[0]->NUMCURRENTSOURCES = CONTACURR;
        // wires
    } while (end_do);

    // ahora slanted

    if (this->swires->n_sw != 0) {
        hay_slanted_wires = true;
    } else {
        hay_slanted_wires = false;
    }

    // SLANTED WIRES
    for (j = 1; j <= this->swires->n_sw; j++) {
        contamedia = contamedia + 1;
        allocate(sgg->Med[contamedia]->SlantedWire[0]);

        sgg->Med[contamedia]->Epr = sgg->Med[1]->Epr;
        sgg->Med[contamedia]->Sigma = sgg->Med[1]->Sigma;
        sgg->Med[contamedia]->Mur = sgg->Med[1]->Mur;
        sgg->Med[contamedia]->SigmaM = sgg->Med[1]->SigmaM;
        sgg->Med[contamedia]->Is.SlantedWire = true;
        sgg->Med[contamedia]->Is.Dielectric = false;
        sgg->Med[contamedia]->SlantedWire[0]->radius = this->swires->SW[j]->RAD;
        if (boundwireradius) {
            if (sgg->Med[contamedia]->SlantedWire[0]->radius > maxwireradius) sgg->Med[contamedia]->SlantedWire[0]->radius = maxwireradius;
        }
        sgg->Med[contamedia]->SlantedWire[0]->R = this->swires->SW[j]->res;
        sgg->Med[contamedia]->SlantedWire[0]->L = this->swires->SW[j]->ind;
        sgg->Med[contamedia]->SlantedWire[0]->C = this->swires->SW[j]->cap;
        sgg->Med[contamedia]->SlantedWire[0]->P_R = this->swires->SW[j]->P_res;
        sgg->Med[contamedia]->SlantedWire[0]->P_L = this->swires->SW[j]->P_ind;
        sgg->Med[contamedia]->SlantedWire[0]->P_C = this->swires->SW[j]->P_cap;

        if (this->swires->SW[j]->disp) {
            allocate(sgg->Med[contamedia]->SlantedWire[0]->disp[0]);
            asignawiredisper(sgg->Med[contamedia]->SlantedWire[0]->disp[0],
                this->swires->SW[j]->dispfile);
        }

        sgg->Med[contamedia]->SlantedWire[0]->LeftEnd = this->swires->SW[j]->LeftEnd;
        sgg->Med[contamedia]->SlantedWire[0]->RightEnd = this->swires->SW[j]->RightEnd;

        sgg->Med[contamedia]->SlantedWire[0]->HasParallel_LeftEnd = false;
        sgg->Med[contamedia]->SlantedWire[0]->Parallel_R_LeftEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Parallel_L_LeftEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Parallel_C_LeftEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->HasParallel_RightEnd = false;
        sgg->Med[contamedia]->SlantedWire[0]->Parallel_R_RightEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Parallel_L_RightEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Parallel_C_RightEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->HasSeries_LeftEnd = false;
        sgg->Med[contamedia]->SlantedWire[0]->Series_R_LeftEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Series_L_LeftEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Series_C_LeftEnd = 2.0e7_RKIND; // en corto 14/2/14
        sgg->Med[contamedia]->SlantedWire[0]->HasSeries_RightEnd = false;
        sgg->Med[contamedia]->SlantedWire[0]->Series_R_RightEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Series_L_RightEnd = 0.0_RKIND;
        sgg->Med[contamedia]->SlantedWire[0]->Series_C_RightEnd = 2.0e7_RKIND; // en corto 14/2/14

        if (this->swires->sw[j]->TL == Parallel_CONS) {
            sgg->Med[contamedia]->SlantedWire[0]->HasParallel_LeftEnd = true;
            sgg->Med[contamedia]->SlantedWire[0]->Parallel_R_LeftEnd = this->swires->sw[j]->R_LeftEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Parallel_L_LeftEnd = this->swires->sw[j]->L_LeftEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Parallel_C_LeftEnd = this->swires->sw[j]->C_LeftEnd;
        } else if (this->swires->sw[j]->TL == SERIES_CONS) {
            sgg->Med[contamedia]->SlantedWire[0]->HasSeries_LeftEnd = true;
            sgg->Med[contamedia]->SlantedWire[0]->Series_R_LeftEnd = this->swires->sw[j]->R_LeftEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Series_L_LeftEnd = this->swires->sw[j]->L_LeftEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Series_C_LeftEnd = this->swires->sw[j]->C_LeftEnd;
        } else if (this->swires->SW[j]->TL == DISPERSIVE_CONS) {
            allocate(sgg->Med[contamedia]->SlantedWire[0]->disp_LeftEnd[0]);
            asignawiredisper(sgg->Med[contamedia]->SlantedWire[0]->disp_LeftEnd[0],
                this->swires->SW[j]->dispfile_LeftEnd);
        }
        if (this->swires->sw[j]->TR == Parallel_CONS) {
            sgg->Med[contamedia]->SlantedWire[0]->HasParallel_RightEnd = true;
            sgg->Med[contamedia]->SlantedWire[0]->Parallel_R_RightEnd = this->swires->sw[j]->R_RightEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Parallel_L_RightEnd = this->swires->sw[j]->L_RightEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Parallel_C_RightEnd = this->swires->sw[j]->C_RightEnd;
        } else if (this->swires->sw[j]->TR == SERIES_CONS) {
            sgg->Med[contamedia]->SlantedWire[0]->HasSeries_RightEnd = true;
            sgg->Med[contamedia]->SlantedWire[0]->Series_R_RightEnd = this->swires->sw[j]->R_RightEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Series_L_RightEnd = this->swires->sw[j]->L_RightEnd;
            sgg->Med[contamedia]->SlantedWire[0]->Series_C_RightEnd = this->swires->sw[j]->C_RightEnd;
        } else if (this->swires->SW[j]->TR == DISPERSIVE_CONS) {
            allocate(sgg->Med[contamedia]->SlantedWire[0]->disp_RightEnd[0]);
            asignawiredisper(sgg->Med[contamedia]->SlantedWire[0]->disp_RightEnd[0],
                this->swires->SW[j]->dispfile_RightEnd);
        }

        sgg->Med[contamedia]->SlantedWire[0]->numNodes = this->swires->SW[j]->n_swc;
        allocate(sgg->Med[contamedia]->SlantedWire[0]->nodes[1..this->swires->SW[j]->n_swc]);
        for (i = 1; i <= this->swires->SW[j]->n_swc; i++) {
            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->VsourceExists = false;
            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->IsourceExists = false;
            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Vsource = nullptr;
            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Isource = nullptr;

            // 2019 clipeo 2 celdas antes
            if (CLIPREGION) {
                if (this->swires->SW[j]->swc[i]->x >= SINPML_Fullsize[iHx]->XE - 2) this->swires->SW[j]->swc[i]->x = SINPML_Fullsize[iHx]->XE - 2;
                if (this->swires->SW[j]->swc[i]->x <= SINPML_Fullsize[iHx]->XI + 2) this->swires->SW[j]->swc[i]->x = SINPML_Fullsize[iHx]->XI + 2;
                if (this->swires->SW[j]->swc[i]->y >= SINPML_Fullsize[iHy]->YE - 2) this->swires->SW[j]->swc[i]->y = SINPML_Fullsize[iHy]->YE - 2;
                if (this->swires->SW[j]->swc[i]->y <= SINPML_Fullsize[iHy]->YI + 2) this->swires->SW[j]->swc[i]->y = SINPML_Fullsize[iHy]->YI + 2;
                if (this->swires->SW[j]->swc[i]->z >= SINPML_Fullsize[iHz]->ZE - 2) this->swires->SW[j]->swc[i]->z = SINPML_Fullsize[iHz]->ZE - 2;
                if (this->swires->SW[j]->swc[i]->z <= SINPML_Fullsize[iHz]->ZI + 2) this->swires->SW[j]->swc[i]->z = SINPML_Fullsize[iHz]->ZI + 2;
            }

            // fin clipeo

            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->index = this->swires->SW[j]->swc[i]->nd;
            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->x = this->swires->SW[j]->swc[i]->x;
            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->y = this->swires->SW[j]->swc[i]->y;
            sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->z = this->swires->SW[j]->swc[i]->z;
            numertag = searchtag(tagtype, this->swires->SW[j]->swc[i]->tag);

            // 2019 para informacion bbox hilos
            // !!!!!!2020 retocado
            if (static_cast<int>(this->swires->SW[j]->swc[i]->x) < bboxwirxi) bboxwirXI = static_cast<int>(this->swires->SW[j]->swc[i]->x);
            if (static_cast<int>(this->swires->SW[j]->swc[i]->x) + 1 > bboxwirxE) bboxwirXE = static_cast<int>(this->swires->SW[j]->swc[i]->x) + 1;
            if (static_cast<int>(this->swires->SW[j]->swc[i]->y) < bboxwirYi) bboxwirYI = static_cast<int>(this->swires->SW[j]->swc[i]->y);
            if (static_cast<int>(this->swires->SW[j]->swc[i]->y) + 1 > bboxwirYE) bboxwirYE = static_cast<int>(this->swires->SW[j]->swc[i]->y) + 1;
            if (static_cast<int>(this->swires->SW[j]->swc[i]->z) < bboxwirZi) bboxwirZI = static_cast<int>(this->swires->SW[j]->swc[i]->z);
            if (static_cast<int>(this->swires->SW[j]->swc[i]->z) + 1 > bboxwirZE) bboxwirZE = static_cast<int>(this->swires->SW[j]->swc[i]->z) + 1;
            // !!!!!

            if (trim(adjustl(this->swires->SW[j]->swc[i]->SRCTYPE)) == F_SOURCE_VOLTAGE) {
                allocate(sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Vsource);
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->VsourceExists = true;
                // sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->VSource->SOFT = true; // fuentes duras sgg 230323. default blandas
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->VSource->Multiplier = this->swires->SW[j]->swc[i]->m;
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->VSource->Resistance = 0.0_RKIND;
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->VSource->fichero.name = trim(adjustl(this->swires->SW[j]->swc[i]->SRCFILE));
            } else if (trim(adjustl(this->swires->SW[j]->swc[i]->SRCTYPE)) == F_SOURCE_CURRENT) {
                allocate(sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Isource);
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->IsourceExists = true;
                // sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Isource->SOFT = true;
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Isource->Multiplier = this->swires->SW[j]->swc[i]->m;
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Isource->Resistance = 0.0_RKIND;
                sgg->Med[contamedia]->SlantedWire[0]->nodes[i]->Isource->fichero.name = trim(adjustl(this->swires->SW[j]->swc[i]->SRCFILE));
            } else if (trim(adjustl(this->swires->SW[j]->swc[i]->SRCTYPE)) != "None") {
                sprintf(buff, "WRONG type of wire source %s", trim(adjustl(this->swires->SW[j]->swc[i]->SRCTYPE)));
                stoponerror(layoutnumber, num_procs, buff);
            }
        }
    }
    // END SLANTED WIRES

    for (j = 1; j <= ubound(conformal_media, 1); j++) {
        addConformalMedia(sgg, media, conformal_media[j], edge_ratios, face_ratios, contamedia, conf_bounding_box, side_to_triangles_maps[j]);
        numertag = searchtag(tagtype, conformal_media[j]->tag);
        CreateConformalPECVolume(layoutnumber, media->sggMtag, tag_numbers, numertag, media->sggMiEx, media->sggMiEy, media->sggMiEz,
            media->sggMiHx, media->sggMiHy, media->sggMiHz, Alloc_iEx_XI,
            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg->Med, sgg->NumMedia, conf_bounding_box, 0);
    }

    contamedia = contamedia + ubound(edge_ratios, 1) + ubound(face_ratios, 1);
    if (findloc(edge_ratios, 0.0, 1) != 0) contamedia = contamedia - 1;
    if (findloc(face_ratios, 0.0, 1) != 0) contamedia = contamedia - 1;

#ifdef CompileWithMTLN
    {
        cable_t* ptr = nullptr;
        for (j = 1; j <= this->mtln->n_sh + this->mtln->n_unsh; j++) {
            ptr = this->mtln->cables[j]->ptr;
            if (dynamic_cast<unshielded_multiwire_t*>(ptr)) {
                contamedia = contamedia + 1;
                allocate(sgg->Med[contamedia]->multiwire[0]);
                sgg->Med[contamedia]->Priority = prior_TW;

                sgg->Med[contamedia]->Epr = sgg->Med[1]->Epr;
                sgg->Med[contamedia]->Sigma = sgg->Med[1]->Sigma;
                sgg->Med[contamedia]->Mur = sgg->Med[1]->Mur;
                sgg->Med[contamedia]->SigmaM = sgg->Med[1]->SigmaM;
                sgg->Med[contamedia]->Is.Multiwire = true;

                isathinwire = false;
                numertag = searchtag(tagtype, this->mtln->cables[j]->ptr->tag);
                for (k = 1; k <= ptr->n_segments; k++) {
                    punto.xi = ptr->segments[k]->x;
                    punto.xe = ptr->segments[k]->x;
                    punto.yi = ptr->segments[k]->y;
                    punto.ye = ptr->segments[k]->y;
                    punto.zi = ptr->segments[k]->z;
                    punto.ze = ptr->segments[k]->z;
                    orientacion = ptr->segments[k]->orientation;
                    CreateLineMM(layoutnumber, media->sggMtag, tag_numbers, numertag, media->sggMiEx, media->sggMiEy, media->sggMiEz,
                        media->sggMiHx, media->sggMiHy, media->sggMiHz, Alloc_iEx_XI,
                        Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                        Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                        Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                        Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                        Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg->Med, sgg->NumMedia, sgg->EShared, BoundingBox, punto, orientacion,
                        contamedia, isathinwire, verbose, numeroasignaciones);
                }
            }
        }
    }
#endif
    // reporta el bounding box

#ifdef CompileWithMPI
    MPI_Barrier(MPI_COMM_WORLD, &ierr);
    MPI_AllReduce(&bboxwirXI, &dummy_bboxwirXI, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, &ierr);
    MPI_AllReduce(&bboxwirYI, &dummy_bboxwirYI, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, &ierr);
    MPI_AllReduce(&bboxwirZI, &dummy_bboxwirzI, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, &ierr);
    MPI_AllReduce(&bboxwirXE, &dummy_bboxwirXE, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &ierr);
    MPI_AllReduce(&bboxwirYE, &dummy_bboxwirYE, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &ierr);
    MPI_AllReduce(&bboxwirZE, &dummy_bboxwirZE, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &ierr);
    MPI_Barrier(MPI_COMM_WORLD, &ierr);
    bboxwirXI = dummy_bboxwirXI;
    bboxwirYI = dummy_bboxwirYI;
    bboxwirzI = dummy_bboxwirzI;
    bboxwirXE = dummy_bboxwirXE;
    bboxwirYE = dummy_bboxwirYE;
    bboxwirZE = dummy_bboxwirZE;
#endif
    sprintf(buff, "pre1_INFO:  Bounding Box for WIREs min_x,min_y,min_z, MAX_x,MAX_y,MAX_z %12d%12d%12d%12d%12d%12d", bboxwirXI, bboxwirYI, bboxwirZI, bboxwirXE, bboxwirYE, bboxwirZE);
    if (((bboxwirXI < (1 << 20)) || (bboxwirYI < (1 << 20)) || (bboxwirZI < (1 << 20)) || (bboxwirXE > -(1 << 20)) || (bboxwirYE > -(1 << 20)) || (bboxwirZE > -(1 << 20))) || VERBOSE) {
        WarnErrReport(buff);
    }
    // FIN WIRES

    if (run_with_dmma) {
        // always at the end since the orientation is found from the PEC one
        // thin Slots
        // the embedding material properties are also if also needed
        tama = this->tSlots->n_tg;
        //
        for (j = 1; j <= tama; j++) {
            //
            //
            tama2 = this->tSlots->Tg[j]->N_tgc;
            for (i = 1; i <= tama2; i++) {
                // del Slot
                direccion = this->tSlots->Tg[j]->TgC[i]->dir;
                i1 = this->tSlots->Tg[j]->TgC[i]->i;
                j1 = this->tSlots->Tg[j]->TgC[i]->j;
                k1 = this->tSlots->Tg[j]->TgC[i]->k;
                ORIX = false;
                ORIY = false;
                ORIZ = false;
                ORIX2 = false;
                ORIY2 = false;
                ORIZ2 = false;
                ORIX3 = false;
                ORIY3 = false;
                ORIZ3 = false;
                ORIX4 = false;
                ORIY4 = false;
                ORIZ4 = false;
                //
                if ((i1 >= BoundingBox.XI) && (i1 < BoundingBox.XE) &&

} else if ((j1 >= BoundingBox.YI) && (j1 < BoundingBox.YE) && (k1 >= BoundingBox.ZI) && (k1 < BoundingBox.ZE)) {
                // encuentra la orientacion del plano PEC que contiene al Slot
                oriX = (direccion == iEy) &&
                       (((media.sggMiHx(i1, j1, k1) == 0) || (sgg.med(media.sggMiHx(i1, j1, k1)).is.pec)) ||
                        (sgg.Med(media.sggMiHx(i1, j1, k1)).Is.ThinSlot));
                
                oriX4 = (direccion == iEz) &&
                        (((media.sggMiHx(i1, j1, k1) == 0) || (sgg.med(media.sggMiHx(i1, j1, k1)).is.pec)) ||
                         (sgg.Med(media.sggMiHx(i1, j1, k1)).Is.ThinSlot));

                oriY = (direccion == iEx) &&
                       (((media.sggMiHy(i1, j1, k1) == 0) || (sgg.med(media.sggMiHy(i1, j1, k1)).is.pec)) ||
                        (sgg.Med(media.sggMiHy(i1, j1, k1)).Is.ThinSlot));

                oriY4 = (direccion == iEz) &&
                        (((media.sggMiHy(i1, j1, k1) == 0) || (sgg.med(media.sggMiHy(i1, j1, k1)).is.pec)) ||
                         (sgg.Med(media.sggMiHy(i1, j1, k1)).Is.ThinSlot));

                oriZ = (direccion == iEx) &&
                       (((media.sggMiHz(i1, j1, k1) == 0) || (sgg.med(media.sggMiHz(i1, j1, k1)).is.pec)) ||
                        (sgg.Med(media.sggMiHz(i1, j1, k1)).Is.ThinSlot));

                oriZ4 = (direccion == iEy) &&
                        (((media.sggMiHz(i1, j1, k1) == 0) || (sgg.med(media.sggMiHz(i1, j1, k1)).is.pec)) ||
                         (sgg.Med(media.sggMiHz(i1, j1, k1)).Is.ThinSlot));

                // encuentra la orientacion del plano PEC que contiene al Slot (considera los vecinos)
                oriX2 = (direccion == iEy) &&
                        (((media.sggMiHx(i1, j1, k1 - 1) == 0) || (sgg.med(media.sggMiHx(i1, j1, k1 - 1)).is.pec)) ||
                         (sgg.Med(media.sggMiHx(i1, j1, k1 - 1)).Is.ThinSlot));

                oriX3 = (direccion == iEz) &&
                        (((media.sggMiHx(i1, j1 - 1, k1) == 0) || (sgg.med(media.sggMiHx(i1, j1 - 1, k1)).is.pec)) ||
                         (sgg.Med(media.sggMiHx(i1, j1 - 1, k1)).Is.ThinSlot));

                oriY2 = (direccion == iEx) &&
                        (((media.sggMiHy(i1, j1, k1 - 1) == 0) || (sgg.med(media.sggMiHy(i1, j1, k1 - 1)).is.pec)) ||
                         (sgg.Med(media.sggMiHy(i1, j1, k1 - 1)).Is.ThinSlot));

                oriY3 = (direccion == iEz) &&
                        (((media.sggMiHy(i1 - 1, j1, k1) == 0) || (sgg.med(media.sggMiHy(i1 - 1, j1, k1)).is.pec)) ||
                         (sgg.Med(media.sggMiHy(i1 - 1, j1, k1)).Is.ThinSlot));

                oriZ2 = (direccion == iEx) &&
                        (((media.sggMiHz(i1, j1 - 1, k1) == 0) || (sgg.med(media.sggMiHz(i1, j1 - 1, k1)).is.pec)) ||
                         (sgg.Med(media.sggMiHz(i1, j1 - 1, k1)).Is.ThinSlot));

                oriZ3 = (direccion == iEy) &&
                        (((media.sggMiHz(i1 - 1, j1, k1) == 0) || (sgg.med(media.sggMiHz(i1 - 1, j1, k1)).is.pec)) ||
                         (sgg.Med(media.sggMiHz(i1 - 1, j1, k1)).Is.ThinSlot));

                if (oriX || oriX4) {
                    orientacion = iEx;
                } else if (oriY || oriY4) {
                    orientacion = iEy;
                } else if (oriZ || oriZ4) {
                    orientacion = iEz;
                } else if (oriX2) {
                    orientacion = iEx;
                    k1 = k1 - 1;
                } else if (oriY2) {
                    orientacion = iEy;
                    k1 = k1 - 1;
                } else if (oriZ2) {
                    orientacion = iEz;
                    j1 = j1 - 1;
                } else if (oriX3) {
                    orientacion = iEx;
                    j1 = j1 - 1;
                } else if (oriY3) {
                    orientacion = iEy;
                    i1 = i1 - 1;
                } else if (oriZ3) {
                    orientacion = iEz;
                    i1 = i1 - 1;
                } else {
                    std::ostringstream buff;
                    buff << "Cannot determine ortientation of the PEC plane with the Slot" << i1 << " " << j1 << " " << k1 << " " << direccion;
                    stoponerror(layoutnumber, num_procs, buff.str());
                }

                this->tSlots.Tg(j).TgC(i).or = orientacion;
                medio2 = -1;
                medio1 = -1;
                switch (std::abs(orientacion)) {
                    case iEx:
                        medio1 = media.sggMiHx(i1, j1, k1);
                        if (i1 > BoundingBox.XI) {
                            medio2 = media.sggMiHx(i1 - 1, j1, k1);
                        } else {
                            medio2 = medio1;
                        }
                        break;
                    case iEy:
                        medio1 = media.sggMiHy(i1, j1, k1);
                        if (j1 > BoundingBox.YI) {
                            medio2 = media.sggMiHy(i1, j1 - 1, k1);
                        } else {
                            medio2 = medio1;
                        }
                        break;
                    case iEz:
                        medio1 = media.sggMiHz(i1, j1, k1);
                        if (k1 > BoundingBox.ZI) {
                            medio2 = media.sggMiHz(i1, j1, k1 - 1);
                        } else {
                            medio2 = medio1;
                        }
                        break;
                }

                if (((sgg.Med(medio1).Is.Dielectric) || (sgg.Med(medio1).Is.Thinslot) || (sgg.Med(medio1).Is.Pec) || (medio1 == 1)) &&
                    ((sgg.Med(medio2).Is.Dielectric) || (sgg.Med(medio2).Is.Thinslot) || (sgg.Med(medio2).Is.Pec) || (medio2 == 1))) {
                    // average adjacent media
                    epr1 = 0.5_RKIND * (sgg.Med(medio1).Epr + sgg.Med(medio2).Epr);
                    mur1 = 0.5_RKIND * (sgg.Med(medio1).Mur + sgg.Med(medio2).Mur);
                } else {
                    std::ostringstream buff;
                    buff << "Media around the Slot are not plain media: " << medio1 << " " << medio2;
                    STOPONERROR(layoutnumber, num_procs, buff.str());
                }
                width = this->tSlots.Tg(j).width;
                if (sgg.NumPlaneWaves == 1) {
                    // assume the incident plane wave if there are planewaves
                    dir(0) = px;
                    dir(1) = py;
                    dir(2) = pz;
                } else {
                    // assume normal incidence
                    switch (std::abs(orientacion)) {
                        case iEx:
                            dir(0) = 1.0_RKIND;
                            dir(1) = 0.0_RKIND;
                            dir(2) = 0.0_RKIND;
                            break;
                        case iEy:
                            dir(0) = 0.0_RKIND;
                            dir(1) = 1.0_RKIND;
                            dir(2) = 0.0_RKIND;
                            break;
                        case iEz:
                            dir(0) = 0.0_RKIND;
                            dir(1) = 0.0_RKIND;
                            dir(2) = 1.0_RKIND;
                            break;
                    }
                }
                dmma_thin_Slot(sgg.dx(i1), sgg.dy(j1), sgg.dz(k1), dir, orientacion, direccion, width, epr1, mur1, EprSlot, MurSlot, eps0, mu0);
                
                indicemedio = contamedia + 1;
                bool iguales = false;
                for (int ii = 1; ii <= contamedia; ++ii) {
                    if (sgg.Med(ii).Is.ThinSlot) {
                        iguales = true;
                        for (int j11 = 1; j11 <= 3; ++j11) {
                            for (int i11 = 1; i11 <= 3; ++i11) {
                                iguales = iguales && (sgg.Med(ii).Anisotropic(1).Epr(i11, j11) == EprSlot(i11, j11)) &&
                                             (sgg.Med(ii).Anisotropic(1).Mur(i11, j11) == MurSlot(i11, j11));
                            }
                        }
                        if (iguales) {
                            indicemedio = ii;
                            break;
                        }
                    }
                }
                if (indicemedio == contamedia + 1) {
                    contamedia = indicemedio;
                    sgg.Med(contamedia).Anisotropic(1).resize(1);
                    sgg.Med(contamedia).Anisotropic(1)[0].Epr = EprSlot;
                    sgg.Med(contamedia).Anisotropic(1)[0].Mur = MurSlot;
                    // lossless
                    sgg.Med(contamedia).Anisotropic(1)[0].Sigma = 0.0_RKIND;
                    sgg.Med(contamedia).Anisotropic(1)[0].SigmaM = 0.0_RKIND;
                    
                    sgg.Med(contamedia).Epr = epr1;
                    sgg.Med(contamedia).Mur = mur1;
                    sgg.Med(contamedia).Sigma = 0.0_RKIND;
                    sgg.Med(contamedia).SigmaM = 0.0_RKIND;
                    sgg.Med(contamedia).Is.Anisotropic = true;
                    // just for signaling
                    sgg.Med(contamedia).Is.ThinSlot = true;
                    sgg.Med(contamedia).Priority = prior_TG;
                }
                
                // record coordinates
                punto.XI = i1;
                punto.XE = i1;
                punto.YI = j1;
                punto.YE = j1;
                punto.ZI = k1;
                punto.ZE = k1;
                numertag = searchtag(tagtype, this->tSlots.Tg(j).TgC(i).tag);
                CreateSurfaceSlotMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE,
                    Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE,
                    Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                    Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE,
                    Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, sgg.HShared, BoundingBox,
                    punto, orientacion, direccion, indicemedio);
            }
        }
    }
    
    // nodalsource
    // precounting
    tama = this->nodsrc.n_nodSrc;
    std::vector<int> contapuntos(tama * (this->nodsrc.n_C2p_max + this->nodsrc.n_C1p_max));
    contapuntos.assign(contapuntos.size(), 0);
    int conta1 = 0;
    for (int i = 1; i <= tama; ++i) {
        int conta2 = 0;
        int tama2 = this->nodsrc.NodalSource(i).n_c1P;
        int tama3 = this->nodsrc.NodalSource(i).n_c2P;
        for (int ii = 1; ii <= tama2; ++ii) {
            punto_s.or = this->nodsrc.NodalSource(i).c1P[ii].or;
            punto_s.XI = this->nodsrc.NodalSource(i).c1P[ii].XI;
            punto_s.XE = this->nodsrc.NodalSource(i).c1P[ii].XE;
            punto_s.YI = this->nodsrc.NodalSource(i).c1P[ii].YI;
            punto_s.YE = this->nodsrc.NodalSource(i).c1P[ii].YE;
            punto_s.ZI = this->nodsrc.NodalSource(i).c1P[ii].ZI;
            punto_s.ZE = this->nodsrc.NodalSource(i).c1P[ii].ZE;
            if ((punto_s.XI <= punto_s.XE) && (punto_s.YI <= punto_s.YE) && (punto_s.ZI <= punto_s.ZE)) {
                conta2++;
            }
        }
        
        for (int ii = 1; ii <= tama3; ++ii) {
            punto_s.or = this->nodsrc.NodalSource(i).c2P[ii].or;
            punto_s.XI = this->nodsrc.NodalSource(i).c2P[ii].XI;
            punto_s.XE = this->nodsrc.NodalSource(i).c2P[ii].XE;
            punto_s.YI = this->nodsrc.NodalSource(i).c2P[ii].YI;
            punto_s.YE = this->nodsrc.NodalSource(i).c2P[ii].YE;
            punto_s.ZI = this->nodsrc.NodalSource(i).c2P[ii].ZI;
            punto_s.ZE = this->nodsrc.NodalSource(i).c2P[ii].ZE;
            if ((punto_s.XI <= punto_s.XE) && (punto_s.YI <= punto_s.YE) && (punto_s.ZI <= punto_s.ZE)) {
                conta2++;
            }
        }
        if (conta2 != 0) {
            conta1++;
            contapuntos[conta1 - 1] = conta2;
        }
    }
    sgg.NumNodalSources = conta1;
    sgg.NodalSource.resize(conta1);
    if (sgg.NumNodalSources != 0) contamedia = contamedia + 1;
    
    conta1 = 0;
    for (int i = 1; i <= tama; ++i) {
        if (contapuntos[i - 1] != 0) {
            conta1++;
            sgg.NodalSource[conta1 - 1].numpuntos = contapuntos[i - 1];
            sgg.NodalSource[conta1 - 1].punto.resize(contapuntos[i - 1]);
            // initialization
            for (int ii = 0; ii < contapuntos[i - 1]; ++ii) {
                sgg.NodalSource[conta1 - 1].punto[ii].or = 0;
                sgg.NodalSource[conta1 - 1].punto[ii].xc = 0.0_RKIND;
                sgg.NodalSource[conta1 - 1].punto[ii].yc = 0.0_RKIND;
                sgg.NodalSource[conta1 - 1].punto[ii].zc = 0.0_RKIND;
                sgg.NodalSource[conta1 - 1].punto[ii].XI = -1;
                sgg.NodalSource[conta1 - 1].punto[ii].XE = -1;
                sgg.NodalSource[conta1 - 1].punto[ii].YI = -1;
                sgg.NodalSource[conta1 - 1].punto[ii].YE = -1;
                sgg.NodalSource[conta1 - 1].punto[ii].ZI = -1;
                sgg.NodalSource[conta1 - 1].punto[ii].ZE = -1;
            }
        }
    }
    
    // asignacion
    conta1 = 0;
    for (int i = 1; i <= tama; ++i) {
        int conta2 = 0;
        if (contapuntos[i - 1] != 0) {
            conta1++;
            sgg.NodalSource[conta1 - 1].fichero.name = trim(adjustl(this->nodsrc.NodalSource(i).nombre));
            sgg.NodalSource[conta1 - 1].isElec = this->nodsrc.NodalSource(i).isElec;
            sgg.NodalSource[conta1 - 1].IsHard = this->nodsrc.NodalSource(i).isHard;
            sgg.NodalSource[conta1 - 1].IsInitialValue = this->nodsrc.NodalSource(i).IsInitialValue;
        }
        
        int tama2 = this->nodsrc.NodalSource(i).n_c1P;
        int tama3 = this->nodsrc.NodalSource(i).n_c2P;
        for (int ii = 1; ii <= tama2; ++ii) {
            // correct bounding box
            punto_s.or = this->nodsrc.NodalSource(i).c1P[ii].or;
            punto_s.xc = this->nodsrc.NodalSource(i).c1P[ii].xc;
            punto_s.yc = this->nodsrc.NodalSource(i).c1P[ii].yc;
            punto_s.zc = this->nodsrc.NodalSource(i).c1P[ii].zc;
            
            punto_s.XI = std::max(this->nodsrc.NodalSource(i).c1P[ii].XI, std::min(BoundingBox.XI, BoundingBox.XE));
            punto_s.YI = std::max(this->nodsrc.NodalSource(i).c1P[ii].YI, std::min(BoundingBox.YI, BoundingBox.YE));
            punto_s.ZI = std::max(this->nodsrc.NodalSource(i).c1P[ii].ZI, std::min(BoundingBox.ZI, BoundingBox.ZE));
            
            punto_s.XE = std::min(this->nodsrc.NodalSource(i).c1P[ii].XE, std::max(BoundingBox.XI, BoundingBox.XE));
            punto_s.YE = std::min(this->nodsrc.NodalSource(i).c1P[ii].YE, std::max(BoundingBox.YI, BoundingBox.YE));
            if ((punto_s.zc != 0) && (this->nodsrc.NodalSource(i).isElec)) { // only in case of Ez
                punto_s.ZE = std::min(this->nodsrc.NodalSource(i).c1P[ii].ZE, std::max(BoundingBox.ZI, BoundingBox.ZE - 1));
            } else {
                punto_s.ZE = std::min(this->nodsrc.NodalSource(i).c1P[ii].ZE, std::max(BoundingBox.ZI, BoundingBox.ZE));
            }
        }
    }

}
        }
        //
        //
        for (k1 = punto_s.ZI; k1 <= punto_s.ZE; k1++) {
            for (j1 = punto_s.YI; j1 <= punto_s.YE; j1++) {
                for (i1 = punto_s.XI; i1 <= punto_s.XE; i1++) {
                    if (punto_s.xc != 0) {
                        //bug OLD 181214 sl_4_20mm_gli.nfde. Fuente nodal electrica embebida en pec y nodal magnetica en pmc se ignoraran sean hard or soft
                        MEDIO = media.sggMiEx(i1, j1, k1);
                        valido = true;
                        //if ( !this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        if (this->nodsrc->NodalSource(i).isElec) {
                            VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PEC);
                        }
                        write(buff, "WARNING: Ex Nodal source on PEC media will be ignored (", i1, j1, k1, ")");
                        if (!VALIDO) WarnErrReport(buff);
                        //
                        // COMENTADO 250816 PQ DA UN ERROR JUISTO CUANDO CAE LA FUENTE EN UN CORTE MPI. HABRIA QU TOCA LA CASUISTICA DE punto_s%ZE = Min (this%nodsrc%NodalSource(i)%c1P(ii)%ZE, Max(BoundingBox%ZI, BoundingBox%ZE  )) PERO NO LO HE QUERIDO HACER
                        //MEDIO = sggmiHx(i1, j1, k1);
                        //valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        //if (!this->nodsrc->NodalSource(i).isElec) {
                        //   VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PMC);
                        //end if
                        //write(buff, "WARNING: Hx Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        //if (!VALIDO) WarnErrReport(buff);
                    }
                    //
                    //
                    if (punto_s.yc != 0) {
                        MEDIO = media.sggMiEy(i1, j1, k1);
                        valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        if (this->nodsrc->NodalSource(i).isElec) {
                            VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PEC);
                        }
                        write(buff, "WARNING: Ey Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        if (!VALIDO) WarnErrReport(buff);
                        //
                        //
                        //MEDIO = sggmiHy(i1, j1, k1);
                        //valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        //if (!this->nodsrc->NodalSource(i).isElec) {
                        //   VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PMC);
                        //end if
                        //write(buff, "WARNING: Hy Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        //if (!VALIDO) WarnErrReport(buff);
                    }
                    //
                    //
                    if (punto_s.zc != 0) {
                        MEDIO = media.sggMiEz(i1, j1, k1);
                        valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        if (this->nodsrc->NodalSource(i).isElec) {
                            VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PEC);
                        }
                        write(buff, "WARNING: Ez Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        if (!VALIDO) WarnErrReport(buff);
                        //
                        //
                        //MEDIO = sggmiHz(i1, j1, k1);
                        //valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        //if (!this->nodsrc->NodalSource(i).isElec) {
                        //   VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PMC);
                        //end if
                        //write(buff, "WARNING: Hz Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        //if (!VALIDO) WarnErrReport(buff);
                    }
                }
            }
        }
        //
        //
        if ((punto_s.XI <= punto_s.XE) && (punto_s.YI <= punto_s.YE) && (punto_s.ZI <= punto_s.ZE)) {
            conta2 = conta2 + 1;
            sgg->NodalSource(conta1).punto(conta2).or = punto_s.or;
            sgg->NodalSource(conta1).punto(conta2).xc = punto_s.xc;
            sgg->NodalSource(conta1).punto(conta2).yc = punto_s.yc;
            sgg->NodalSource(conta1).punto(conta2).zc = punto_s.zc;
            sgg->NodalSource(conta1).punto(conta2).XI = punto_s.XI;
            sgg->NodalSource(conta1).punto(conta2).XE = punto_s.XE;
            sgg->NodalSource(conta1).punto(conta2).YI = punto_s.YI;
            sgg->NodalSource(conta1).punto(conta2).YE = punto_s.YE;
            sgg->NodalSource(conta1).punto(conta2).ZI = punto_s.ZI;
            sgg->NodalSource(conta1).punto(conta2).ZE = punto_s.ZE;
            //PARA ACOMODAR LAS NODAL SOURCE COMO MEDIOS LINE Y PODER VISUALIZAR SONDAS 010824
            sgg->Med(contamedia).Is.Dielectric = true;
            sgg->Med(contamedia).Is.LINE = true;
            sgg->Med(contamedia).Priority = prior_IL;
            sgg->Med(contamedia).Epr = 1.0;
            sgg->Med(contamedia).Sigma = 0.;
            sgg->Med(contamedia).Mur = 1.0;
            sgg->Med(contamedia).SigmaM = 0.;
            punto.XI = punto_s.XI;
            punto.XE = punto_s.XE;
            punto.YI = punto_s.YI;
            punto.YE = punto_s.YE;
            punto.ZI = punto_s.ZI;
            punto.ZE = punto_s.ZE;
            orientacion = punto_s.or;
            isathinwire = false;
            numertag = 37;
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg->Med, sgg->NumMedia, sgg->EShared, BoundingBox, punto, orientacion,
                contamedia, isathinwire, verbose, numeroasignaciones);
        }
        sgg->NodalSource(conta1).numpuntos = conta2; //update with the correct value
    }
    //
    //
    for (ii = 1; ii <= tama3; ii++) {
        punto_s.or = this->nodsrc->NodalSource(i).c2P(ii).or;
        punto_s.XI = this->nodsrc->NodalSource(i).c2P(ii).XI;
        punto_s.XE = this->nodsrc->NodalSource(i).c2P(ii).XE;
        punto_s.YI = this->nodsrc->NodalSource(i).c2P(ii).YI;
        punto_s.YE = this->nodsrc->NodalSource(i).c2P(ii).YE;
        punto_s.ZI = this->nodsrc->NodalSource(i).c2P(ii).ZI;
        punto_s.ZE = this->nodsrc->NodalSource(i).c2P(ii).ZE;
        punto_s.xc = this->nodsrc->NodalSource(i).c2P(ii).xc;
        punto_s.yc = this->nodsrc->NodalSource(i).c2P(ii).yc;
        punto_s.zc = this->nodsrc->NodalSource(i).c2P(ii).zc;
        //correct bounding box
        punto_s.XI = Max(this->nodsrc->NodalSource(i).c2p(ii).XI, Min(BoundingBox.XI, BoundingBox.XE));
        punto_s.YI = Max(this->nodsrc->NodalSource(i).c2p(ii).YI, Min(BoundingBox.YI, BoundingBox.YE));
        punto_s.ZI = Max(this->nodsrc->NodalSource(i).c2p(ii).ZI, Min(BoundingBox.ZI, BoundingBox.ZE));
        //
        punto_s.XE = Min(this->nodsrc->NodalSource(i).c2p(ii).XE, Max(BoundingBox.XI, BoundingBox.XE));
        punto_s.YE = Min(this->nodsrc->NodalSource(i).c2p(ii).YE, Max(BoundingBox.YI, BoundingBox.YE));
        if ((punto_s.zc != 0) && (this->nodsrc->NodalSource(i).isElec)) { //only in case of Ez
            punto_s.ZE = Min(this->nodsrc->NodalSource(i).c2p(ii).ZE, Max(BoundingBox.ZI, BoundingBox.ZE - 1));
        } else {
            punto_s.ZE = Min(this->nodsrc->NodalSource(i).c2p(ii).ZE, Max(BoundingBox.ZI, BoundingBox.ZE));
        }
        //
        punto_s.or = this->nodsrc->NodalSource(i).c2p(ii).or;
        punto_s.xc = this->nodsrc->NodalSource(i).c2p(ii).xc;
        punto_s.yc = this->nodsrc->NodalSource(i).c2p(ii).yc;
        punto_s.zc = this->nodsrc->NodalSource(i).c2p(ii).zc;
        //
        //
        for (k1 = punto_s.ZI; k1 <= punto_s.ZE; k1++) {
            for (j1 = punto_s.YI; j1 <= punto_s.YE; j1++) {
                for (i1 = punto_s.XI; i1 <= punto_s.XE; i1++) {
                    if (punto_s.xc != 0) {
                        MEDIO = media.sggMiEx(i1, j1, k1);
                        valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        if (this->nodsrc->NodalSource(i).isElec) {
                            VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PEC);
                        }
                        write(buff, "WARNING: Ex Nodal source on PEC media will be ignored (", i1, j1, k1, ")");
                        if (!VALIDO) WarnErrReport(buff);
                        //
                        //
                        //
                        //MEDIO = sggmiHx(i1, j1, k1);
                        //valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        //if (!this->nodsrc->NodalSource(i).isElec) {
                        //   VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PMC);
                        //end if
                        //write(buff, "WARNING: Hx Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        //if (!VALIDO) WarnErrReport(buff);
                    }
                    //
                    //
                    if (punto_s.yc != 0) {
                        MEDIO = media.sggMiEy(i1, j1, k1);
                        valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        if (this->nodsrc->NodalSource(i).isElec) {
                            VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PEC);
                        }
                        write(buff, "WARNING: Ey Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        if (!VALIDO) WarnErrReport(buff);
                        //
                        //
                        //MEDIO = sggmiHy(i1, j1, k1);
                        //valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        //if (!this->nodsrc->NodalSource(i).isElec) {
                        //   VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PMC);
                        //end if
                        //write(buff, "WARNING: Hy Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        //if (!VALIDO) WarnErrReport(buff);
                    }
                    //
                    //
                    if (punto_s.zc != 0) {
                        MEDIO = media.sggMiEz(i1, j1, k1);
                        valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        if (this->nodsrc->NodalSource(i).isElec) {
                            VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PEC);
                        }
                        write(buff, "WARNING: Ez Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        if (!VALIDO) WarnErrReport(buff);
                        //
                        //
                        //MEDIO = sggmiHz(i1, j1, k1);
                        //valido = true;
                        //if (!this->nodsrc->NodalSource(i).isHard) {
                        //  VALIDO = (sgg->Med(MEDIO).Is.Dielectric) || (sgg->Med(MEDIO).Is.EDispersive) || &
                        // & (sgg->Med(MEDIO).Is.MDispersive)
                        //ELSE
                        //  VALIDO = true;
                        //end if
                        //if (!this->nodsrc->NodalSource(i).isElec) {
                        //   VALIDO = VALIDO && (!sgg->Med(MEDIO).Is.PMC);
                        //end if
                        //write(buff, "WARNING: Hz Nodal source on PMC media will be ignored (", i1, j1, k1, ")");
                        //if (!VALIDO) WarnErrReport(buff);
                    }
                }
            }
        }
        //
        //
        if ((punto_s.XI <= punto_s.XE) && (punto_s.YI <= punto_s.YE) && (punto_s.ZI <= punto_s.ZE)) {
            conta2 = conta2 + 1;
            sgg->NodalSource(conta1).punto(conta2).or = punto_s.or;
            sgg->NodalSource(conta1).punto(conta2).xc = punto_s.xc;
            sgg->NodalSource(conta1).punto(conta2).yc = punto_s.yc;
            sgg->NodalSource(conta1).punto(conta2).zc = punto_s.zc;
            sgg->NodalSource(conta1).punto(conta2).XI = punto_s.XI;
            sgg->NodalSource(conta1).punto(conta2).XE = punto_s.XE;
            sgg->NodalSource(conta1).punto(conta2).YI = punto_s.YI;
            sgg->NodalSource(conta1).punto(conta2).YE = punto_s.YE;
            sgg->NodalSource(conta1).punto(conta2).ZI = punto_s.ZI;
            sgg->NodalSource(conta1).punto(conta2).ZE = punto_s.ZE;
            //PARA ACOMODAR LAS NODAL SOURCE COMO MEDIOS LINE Y PODER VISUALIZAR SONDAS 010824
            sgg->Med(contamedia).Is.Dielectric = true;
            sgg->Med(contamedia).Is.LINE = true;
            sgg->Med(contamedia).Priority = prior_IL;
            sgg->Med(contamedia).Epr = 1.0;
            sgg->Med(contamedia).Sigma = 0.;
            sgg->Med(contamedia).Mur = 1.0;
            sgg->Med(contamedia).SigmaM = 0.;
            punto.XI = punto_s.XI;
            punto.XE = punto_s.XE;
            punto.YI = punto_s.YI;
            punto.YE = punto_s.YE;
            punto.ZI = punto_s.ZI;
            punto.ZE = punto_s.ZE;
            orientacion = punto_s.or;
            isathinwire = false;
            numertag = 37;
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg->Med, sgg->NumMedia, sgg->EShared, BoundingBox, punto, orientacion,
                contamedia, isathinwire, verbose, numeroasignaciones);
        }
        sgg->NodalSource(conta1).numpuntos = conta2; //update with the correct value
    }
}
//
if (contapuntos.is_allocated()) contapuntos.deallocate();

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!PROBES
//!!!!!!!!!!!!!!!!!!!!PROBES
//!!!!!!!!!!!!!!!!!!!!PROBES
//!!!!!!!!!!!!!!!!!!!!PROBES
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRIMERO LAS CUENTO
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MasSondas     fields
//
tamaScrPrb = this->VolPrb.length;
if (createmapvtk) tamaScrPrb = tamaScrPrb + 1;
tamaScrPrb = tamaScrPrb * 3; //210618 allocateo el trip

// le por si las volumicas son de tipo tifr mezcladas
tamaSonda = this->Sonda.length;
tamaoldSONDA = this->oldSONDA.n_probes;
tamaBloquePrb = this->BloquePRB.N_BP;
//probes totales
sgg.NumberRequest = (tamaSonda) + (tamaoldSONDA) + (tamaBloquePrb) + (tamaScrPrb);
allocate(sgg.observation, 1, sgg.NumberRequest);
//inicializacion
for (int idx = 1; idx <= sgg.NumberRequest; ++idx) {
    sgg.observation[idx].nP = -1;
    sgg.observation[idx].InitialTime = -1;
    sgg.observation[idx].FinalTime = -1;
    sgg.observation[idx].TimeStep = -1;
    sgg.observation[idx].InitialFreq = -1.0_RKIND;
    sgg.observation[idx].FinalFreq = -1.0_RKIND;
    sgg.observation[idx].FreqStep = -1.0_RKIND;
    sgg.observation[idx].outputrequest = ' ';
    sgg.observation[idx].FileNormalize = ' ';
    sgg.observation[idx].FreqDomain = false;
    sgg.observation[idx].TimeDomain = false;
    sgg.observation[idx].Saveall = false;
    sgg.observation[idx].TRANSFER = false;
    sgg.observation[idx].Volumic = false;
}

//
//ahora las cuento por bloques
//
for (int i = 1; i <= tamaSonda; ++i) {
    ii = i;
    sgg.observation[ii].nP = 0;
    tama2 = (this->Sonda.collection[i].len_cor);
    for (int j = 1; tama2; ++j) {
        tipotemp = this->Sonda.collection[i].cordinates[j].or;
        punto.XI = this->Sonda.collection[i].cordinates[j].XI;
        punto.YI = this->Sonda.collection[i].cordinates[j].YI;
        punto.ZI = this->Sonda.collection[i].cordinates[j].ZI;
        if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) && 
            (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) && 
            (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE) && 
            ((tipotemp == NP_COR_EX) || (tipotemp == NP_COR_EY) || (tipotemp == NP_COR_EZ) || 
             (tipotemp == NP_COR_HX) || (tipotemp == NP_COR_HY) || (tipotemp == NP_COR_HZ))) {
            sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
        } else if (tipotemp == NP_COR_WIRECURRENT || tipotemp == NP_COR_CHARGE) {
            nodo_cazado = false;
            for (int j1 = 1; j1 <= this->twires.n_tw; ++j1) {
                for (int i1 = 1; i1 <= this->twires.TW[j1].N_TWC; ++i1) {
                    //nodo cazado
                    if (this->twires.TW[j1].TWC[i1].nd == this->Sonda.collection[i].cordinates[j].XI) {
                        punto.XI = this->twires.TW[j1].TWC[i1].i;
                        punto.YI = this->twires.TW[j1].TWC[i1].j;
                        punto.ZI = this->twires.TW[j1].TWC[i1].k;
                        nodo_cazado = true;
                        //
                        if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) && 
                            (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) && 
                            (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                            //
                            switch (this->twires.TW[j1].TWC[i1].D) {
                                case iEx:
                                case iEy:
                                case iEz:
                                    sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                                    break;
                            }
                            goto loop_busqueda1_exit;
                        }
                    }
                }
            }
            loop_busqueda1_exit:;
            //si no lo ha cazado... probamos SLANTED
            if (!nodo_cazado) {
                for (int j1 = 1; j1 <= this->swires.n_sw; ++j1) {
                    for (int i1 = 1; i1 <= this->swires.SW[j1].N_SWC; ++i1) {
                        //nodo cazado
                        if (this->swires.SW[j1].SWC[i1].nd == this->Sonda.collection[i].cordinates[j].XI) {
                            punto.XI = floor(this->swires.SW[j1].SWC[i1].x);
                            punto.YI = floor(this->swires.SW[j1].SWC[i1].y);
                            punto.ZI = floor(this->swires.SW[j1].SWC[i1].z);
                            nodo_cazado = true;
                            //
                            if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) && 
                                (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) && 
                                (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                                sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                                //
                                goto loop_busqueda2_exit;
                            }
                        }
                    }
                }
                loop_busqueda2_exit:;
            }
            //si no lo ha cazado
            if (!nodo_cazado) {
                sprintf(buff, "Current probe not found in WIRE segment %i", this->Sonda.collection[i].cordinates[j].XI);
                StopOnError(layoutnumber, num_procs, buff);
            }
        } else if (this->Sonda.collection[i].cordinates[j].or == NP_COR_DDP) {
            if (run_with_dmma) {
                nodo_cazado = false;
                for (int j1 = 1; j1 <= this->tSlots.n_tg; ++j1) {
                    for (int i1 = 1; i1 <= this->tSlots.Tg[j1].N_tgc; ++i1) {
                        //nodo cazado
                        if (this->tSlots.Tg[j1].TgC[i1].node == this->Sonda.collection[i].cordinates[j].XI) {
                            punto.XI = this->tSlots.Tg[j1].TgC[i1].i;
                            punto.YI = this->tSlots.Tg[j1].TgC[i1].j;
                            punto.ZI = this->tSlots.Tg[j1].TgC[i1].k;
                            nodo_cazado = true;
                            if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) && 
                                (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) && 
                                (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                                sgg.observation[i].nP = sgg.observation[i].nP + 1;
                                goto do_loop_busquedatg1_exit;
                            }
                        }
                    }
                }
                do_loop_busquedatg1_exit:;
                //si no lo ha cazado
                if (!nodo_cazado) {
                    sprintf(buff, "Voltage probe not found %i", this->Sonda.collection[i].cordinates[j].XI);
                    StopOnError(layoutnumber, num_procs, buff);
                }
            } else {
                sprintf(buff, "ERROR: Voltage probe in gaps only available under -dmma flag");
                StopOnError(layoutnumber, num_procs, buff);
            }
        } else if (abs(tipotemp) == NP_COR_LINE) {
            sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
        }
    }
}

//Sondas
//
for (int i = 1; i <= tamaoldSONDA; ++i) {
    //acumulador
    ii = i + tamaSonda;
    //far fields (no es time domain pero una forma especial de ellos)
    tama2 = (this->oldSONDA.probes[i].n_FarField);
    if (tama2 > 1) {
        buff = "Only one Far Field probe allowed";
        STOPONERROR(layoutnumber, num_procs, buff);
    }
    if (tama2 > 0) sgg.observation[ii].nP = 1; //un punto para todo el farfield (es simbolico)
    //electric FIELDS
    tama2 = (this->oldSONDA.probes[i].n_Electric);
    for (int j = 1; j <= tama2; ++j) {
        tama3 = (this->oldSONDA.probes[i].Electric[j].probe.n_cord);
        tama4 = (this->oldSONDA.probes[i].Electric[j].probe.n_cord);
        tama5 = (this->oldSONDA.probes[i].Electric[j].probe.n_cord);
        tama6 = (this->oldSONDA.probes[i].Electric[j].probe.n_cord);
        buff = "TAMANIOS DE PROBES EH RAROS";
        if ((tama3 != tama4) || (tama3 != tama5) || (tama3 != tama6)) {
            STOPONERROR(layoutnumber, num_procs, buff);
        }
        for (int k = 1; k <= tama3; ++k) {
            punto.XI = this->oldSONDA.probes[i].Electric[j].probe.i[k];
            punto.YI = this->oldSONDA.probes[i].Electric[j].probe.j[k];
            punto.ZI = this->oldSONDA.probes[i].Electric[j].probe.k[k];
            if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) && 
                (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) && 
                (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                sgg.observation[ii].nP = sgg.observation[ii].nP + 3;
            }
        }
    }
    //MAGNETIC FIELDS
    tama2 = (this->oldSONDA.probes[i].n_Magnetic);
    for (int j = 1; j <= tama2; ++j) {
        tama3 = (this->oldSONDA.probes[i].Magnetic[j].probe.n_cord);
        tama4 = (this->oldSONDA.probes[i].Magnetic[j].probe.n_cord);
        tama5 = (this->oldSONDA.probes[i].Magnetic[j].probe.n_cord);
        tama6 = (this->oldSONDA.probes[i].Magnetic[j].probe.n_cord);
        buff = "pre1_ERROR: TAMANIOS DE PROBES EH RAROS";
        if ((tama3 != tama4) || (tama3 != tama5) || (tama3 != tama6)) {
            STOPONERROR(layoutnumber, num_procs, buff);
        }
        for (int k = 1; k <= tama3; ++k) {
            punto.XI = this->oldSONDA.probes[i].Magnetic[j].probe.i[k];
            punto.YI = this->oldSONDA.probes[i].Magnetic[j].probe.j[k];
            punto.ZI = this->oldSONDA.probes[i].Magnetic[j].probe.k[k];
            if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) && 
                (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) && 
                (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                sgg.observation[ii].nP = sgg.observation[ii].nP + 3;
            }
        }
    }
}
//Bloque Current Probes
//
for (int i = 1; i <= tamaBloquePrb; ++i) {
    ii = i + tamaSonda + tamaoldSONDA;
    sgg.observation[ii].nP = 0;
    punto.XI = this->BloquePRB.BP[i].i1;
    punto.YI = this->BloquePRB.BP[i].j1;
    punto.ZI = this->BloquePRB.BP[i].k1;
    punto.XE = this->BloquePRB.BP[i].I2;
    punto.YE = this->BloquePRB.BP[i].J2;
    punto.ZE = this->BloquePRB.BP[i].K2;
    !!!
    if (((punto.XI >= BoundingBox.XI) || (punto.XI <= BoundingBox.XE)) && 
        ((punto.YI >= BoundingBox.YI) || (punto.YI <= BoundingBox.YE)) && 
        ((punto.ZI >= BoundingBox.ZI) || (punto.ZI <= BoundingBox.ZE)) && 
        ((punto.XE >= BoundingBox.XI) || (punto.XE <= BoundingBox.XE)) && 
        ((punto.YE >= BoundingBox.YI) || (punto.YE <= BoundingBox.YE)) && 
        ((punto.ZE >= BoundingBox.ZI) || (punto.ZE <= BoundingBox.ZE))) {
        switch (this->BloquePRB.BP[i].NML) {
            case iEx:
                for (int k = this->BloquePRB.BP[i].i1; k <= this->BloquePRB.BP[i].I2; k += this->BloquePRB.BP[i].skip) {
                    sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                }
                break;
            case iEy:
                for (int k = this->BloquePRB.BP[i].j1; k <= this->BloquePRB.BP[i].J2; k += this->BloquePRB.BP[i].skip) {
                    sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                }
                break;
            case iEz:
                for (int k = this->BloquePRB.BP[i].k1; k <= this->BloquePRB.BP[i].K2; k += this->BloquePRB.BP[i].skip) {
                    sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                }
                break;
        }
    }
    //DEL TAMA DEL Bloque CURRENT PROBES
}

//Volumic probes (similar to MasSondas PERO CON PUNTOS FINALES COMO LAS Bloque PROBES)
//ahora las cuento por bloques
for (int i = 1; i <= tamaScrPrb / 3; ++i) { // !!!210618 En realidad hay un tercio
    ii = i + tamaSonda + tamaoldSONDA + tamaBloquePrb;
    sgg.observation[ii].nP = 0;
    // crea una sonda vtk vacion con el instante incial 0 a efectos de mapa
    if (createmapvtk && (i == tamaScrPrb / 3)) { // !!!210618 En realidad hay un tercio
        tama2 = 1;
        for (int j = 1; j <= tama2; ++j) {
            tipotemp = mapvtk;
            punto.XI = SINPML_fullsize[iHx].XI; // !!! +1   !ojo si se cambia aqui tambien mas abajo !le quito 1 para que con condiciones PEC no las pinte !habria que manejar el dibujo de las condiciones aparte
            punto.YI = SINPML_fullsize[iHy].YI; // !!! +1
            punto.ZI = SINPML_fullsize[iHz].ZI; // !!! +1
            punto.XE = SINPML_fullsize[iHx].XE; // !!! -1
            punto.YE = SINPML_fullsize[iHy].YE; // !!! -1
            punto.ZE = SINPML_fullsize[iHz].ZE; // !!! -1
            //!!!               print *,layoutnumber,punto%XI,punto%YI,punto%ZI,punto%XE,punto%YE,punto%ZE
            if (((punto.XI >= BoundingBox.XI) || (punto.XI <= BoundingBox.XE)) && 
                ((punto.YI >= BoundingBox.YI) || (punto.YI <= BoundingBox.YE)) && 
                ((punto.ZI >= BoundingBox.ZI) || (punto.ZI <= BoundingBox.ZE)) && 
                ((punto.XE >= BoundingBox.XI) || (punto.XE <= BoundingBox.XE)) && 
                ((punto.YE >= BoundingBox.YI) || (punto.YE <= BoundingBox.YE)) && 
                ((punto.ZE >= BoundingBox.ZI) || (punto.ZE <= BoundingBox.ZE))) {
                //!!!                   print *,'----Dentro->',layoutnumber,punto%XI,punto%YI,punto%ZI,punto%XE,punto%YE,punto%ZE
                sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
            }
        }
    } else {
        tama2 = (this->VolPrb.collection[i].len_cor);
        for (int j = 1; j <= tama2; ++j) {
            tipotemp = this->VolPrb.collection[i].cordinates[j].or;
            punto.XI = this->VolPrb.collection[i].cordinates[j].XI;
            punto.YI = this->VolPrb.collection[i].cordinates[j].YI;
            punto.ZI = this->VolPrb.collection[i].cordinates[j].ZI;
            punto.XE = this->VolPrb.collection[i].cordinates[j].XE;
            punto.YE = this->VolPrb.collection[i].cordinates[j].YE;
            punto.ZE = this->VolPrb.collection[i].cordinates[j].ZE;

            if (((punto.XI >= BoundingBox.XI) || (punto.XI <= BoundingBox.XE)) && 
                ((punto.YI >= BoundingBox.YI) || (punto.YI <= BoundingBox.YE)) && 
                ((punto.ZI >= BoundingBox.ZI) || (punto.ZI <= BoundingBox.ZE)) && 
                ((punto.XE >= BoundingBox.XI) || (punto.XE <= BoundingBox.XE)) && 
                ((punto.YE >= BoundingBox.YI) || (punto.YE <= BoundingBox.YE)) && 
                ((punto.ZE >= BoundingBox.ZI) || (punto.ZE <= BoundingBox.ZE))) {

                sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
            }
        }
    }
}
// si se lanza con -mapvtk se crea una slice probe para ver la estructura

//Ahora creo los puntos de observacion
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sondas = 0;
for (int ii = 1; ii <= sgg.NumberRequest; ++ii) {
    sondas = sondas + sgg.observation[ii].nP;
}
//
if (sondas > Maxprobes) {
    sprintf(buff, "Too many probes= %i Either  or reduce the number of probes below %i", sondas, Maxprobes);
    STOPONERROR(layoutnumber, num_procs, buff);
}
if (sondas > 1024) {
    sprintf(buff, "Number of probes %i is over 1024. UNIX limit is typically 1024 Check! It will crash if ulimit -n is not set > %i", sondas, sondas);
    WarnErrReport(buff);
}

//luego chequeo las memoria de las sondas si se van de memoria en observation.f90
//      if (sondas*BuffObse*4 > MaxMemoryProbes) then
//        write(buff,*) 'Too much memory for the probes= ', sondas * BuffObse * 4, 'Probes= ', sondas,   &
//         &     'Either reduce the number o&
//       &f probes or recompile decreasing BuffObse ', BuffObse, 'or increasing ', MaxMemoryProbes
//        call STOPONERROR(layoutnumber,num_procs,buff)
//      end if
//
if (sgg.NumberRequest != 0) {
    //alocateo
    for (int ii = 1; ii <= sgg.NumberRequest; ++ii) {
        allocate(sgg.observation[ii].P, 1, sgg.observation[ii].nP);
        for (int idx = 1; idx <= sgg.observation[ii].nP; ++idx) {
            sgg.observation[ii].P[idx].what = nothing; //bug peligroso 2012
        }
        sgg.observation[ii].TimeDomain = false;
        sgg.observation[ii].FreqDomain = false;
        sgg.observation[ii].TRANSFER = false;
        sgg.observation[ii].Volumic = false;
        sgg.observation[ii].FileNormalize = ' ';
        //trancos
        for (int idx = 1; idx <= sgg.observation[ii].nP; ++idx) {
            sgg.observation[ii].P[idx].Xtrancos = 1; //default
            sgg.observation[ii].P[idx].Ytrancos = 1; //default
            sgg.observation[ii].P[idx].Ztrancos = 1; //default
        }
        //fin trancos
        //al final debe quedar igual
        //lo reseteo porque lo reutilizo de contador
        sgg.observation[ii].nP = 0; //reset it
        //
    }
    for (int i = 1; i <= tamaSonda; ++i) {
        ii = i;
        tama2 = (this->Sonda.collection[i].len_cor);
        //
        switch (this->Sonda.collection[i].type2) {
            case NP_T2_time:
                sgg.observation[ii].TimeDomain = true;
                sgg.observation[ii].FreqDomain = false;
                sgg.observation[ii].TRANSFER = false;
                break;
            case NP_T2_FREQ:
                //I will output everything in time and transform it later
                sgg.observation[ii].TimeDomain = true;
                sgg.observation[ii].FreqDomain = true;
                sgg.observation[ii].TRANSFER = false;
                //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                break;
            case NP_T2_TRANSFER:
                sgg.observation[ii].TimeDomain = true;
                sgg.observation[ii].FreqDomain = true;
                sgg.observation[ii].TRANSFER = true;
                buff = "Transfer function only in Frequency Domain";
                //!!           call STOPONERROR(layoutnumber,num_procs,buff)
                break;
            case NP_T2_TIMEFREQ:
                sgg.observation[ii].TimeDomain = true;
                sgg.observation[ii].FreqDomain = true;
                sgg.observation[ii].TRANSFER = false;
                //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                break;
            case NP_T2_TIMETRANSF:
                sgg.observation[ii].TimeDomain = true;
                sgg.observation[ii].FreqDomain = true;
                sgg.observation[ii].TRANSFER = true;
                //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                break;
            case NP_T2_FREQTRANSF:
                sgg.observation[ii].TimeDomain = true;
                sgg.observation[ii].FreqDomain = true;
                sgg.observation[ii].TRANSFER = true;
                //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                break;
            case NP_T2_TIMEFRECTRANSF:
                sgg.observation[ii].TimeDomain = true;
                sgg.observation[ii].FreqDomain = true;
                sgg.observation[ii].TRANSFER = true;
                break;
        }
        //repair info
        //
        if (sgg.observation[ii].FreqDomain) {
            //save everything for later transform
            sgg.observation[ii].Saveall = true;
        }
        if (this->Sonda.collection[i].type1 == NP_T1_AMBOS) {
            sgg.observation[ii].Saveall = true;
            sgg.observation[ii].TimeDomain = true;
        } else {
            // continue
        }
        //
        sprintf(probenumber, "%07i", ii);
        //        sgg%observation(ii)%outputrequest=trim(adjustl(probenumber))//'_'// &
        //                                          trim(adjustl(this%Sonda%collection(i)%output
    }
}

request));
            sgg.observation[ii].outputrequest = trim(adjustl(this.Sonda.collection[i].outputrequest));
            sgg.observation[ii].InitialTime = this.Sonda.collection[i].tstart;
            sgg.observation[ii].FinalTime = this.Sonda.collection[i].tstop;
            sgg.observation[ii].TimeStep = this.Sonda.collection[i].tstep;
            if ((sgg.observation[ii].FinalTime < TINY(1.0_RKIND)) || (sgg.observation[ii].TimeStep < TINY(1.0_RKIND))) 
                sgg.observation[ii].Saveall = true;
            sgg.observation[ii].InitialFreq = this.Sonda.collection[i].fstart;
            sgg.observation[ii].FinalFreq = this.Sonda.collection[i].fstop;
            sgg.observation[ii].FreqStep = this.Sonda.collection[i].fstep;
            sgg.observation[ii].FileNormalize = trim(adjustl(this.Sonda.collection[i].filename));
            //
            if ((sgg.observation[ii].InitialFreq < 0.0) ||
                (sgg.observation[ii].FinalFreq <= 1e-9) ||
                (sgg.observation[ii].FreqStep <= 1e-9)) {
                write(buff, *) "ERROR: Some incorrect frequency domain parameters (initial,final,step) ", sgg.observation[ii].InitialFreq, sgg.observation[ii].FinalFreq, sgg.observation[ii].FreqStep;
                if (sgg.observation[ii].FreqDomain) STOPONERROR(layoutnumber, num_procs, buff);
            }
            //
            for (j = 1; j <= tama2; ++j) {
                punto.XI = this.Sonda.collection[i].cordinates[j].XI;
                punto.YI = this.Sonda.collection[i].cordinates[j].YI;
                punto.ZI = this.Sonda.collection[i].cordinates[j].ZI;
                if (this.Sonda.collection[i].cordinates[j].or == NP_COR_EX) {
                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                        sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                        sgg.observation[ii].P[sgg.observation[ii].nP].XI = punto.XI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].YI = punto.YI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].ZI = punto.ZI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].What = iEx;
                    }
                } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_EY) {
                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                        sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                        sgg.observation[ii].P[sgg.observation[ii].nP].XI = punto.XI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].YI = punto.YI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].ZI = punto.ZI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].What = iEy;
                    }
                } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_EZ) {
                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                        sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                        sgg.observation[ii].P[sgg.observation[ii].nP].XI = punto.XI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].YI = punto.YI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].ZI = punto.ZI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].What = iEz;
                    }
                } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_HX) {
                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                        sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                        sgg.observation[ii].P[sgg.observation[ii].nP].XI = punto.XI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].YI = punto.YI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].ZI = punto.ZI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].What = iHx;
                    }
                } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_HY) {
                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                        sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                        sgg.observation[ii].P[sgg.observation[ii].nP].XI = punto.XI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].YI = punto.YI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].ZI = punto.ZI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].What = iHy;
                    }
                } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_HZ) {
                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                        sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                        sgg.observation[ii].P[sgg.observation[ii].nP].XI = punto.XI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].YI = punto.YI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].ZI = punto.ZI;
                        sgg.observation[ii].P[sgg.observation[ii].nP].What = iHz;
                    }
                } else if ((this.Sonda.collection[i].cordinates[j].or == NP_COR_WIRECURRENT) ||
                           (this.Sonda.collection[i].cordinates[j].or == NP_COR_CHARGE)) {
                    nodo_cazado = false;
                    do {
                        for (j1 = 1; j1 <= this.twires.n_tw; ++j1) {
                            for (i1 = 1; i1 <= this.twires.TW[j1].N_TWC; ++i1) {
                                //nodo cazado
                                if (this.twires.TW[j1].TWC[i1].nd == this.Sonda.collection[i].cordinates[j].XI) {
                                    nodo_cazado = true;
                                    punto.XI = this.twires.TW[j1].TWC[i1].i;
                                    punto.YI = this.twires.TW[j1].TWC[i1].j;
                                    punto.ZI = this.twires.TW[j1].TWC[i1].k;
                                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                                        switch (this.twires.TW[j1].TWC[i1].D) {
                                        case iEx:
                                            sgg.observation[i].nP = sgg.observation[i].nP + 1;
                                            sgg.observation[i].P[sgg.observation[i].nP].node = this.twires.TW[j1].TWC[i1].nd;
                                            sgg.observation[i].P[sgg.observation[i].nP].XI = punto.XI;
                                            sgg.observation[i].P[sgg.observation[i].nP].YI = punto.YI;
                                            sgg.observation[i].P[sgg.observation[i].nP].ZI = punto.ZI;
                                            if (this.Sonda.collection[i].cordinates[j].or == NP_COR_WIRECURRENT) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iJx;
                                            } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_CHARGE) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iQx;
                                            }
                                            //se nota con un indice distinto
                                            //
                                            break;
                                        case iEy:
                                            sgg.observation[i].nP = sgg.observation[i].nP + 1;
                                            sgg.observation[i].P[sgg.observation[i].nP].node = this.twires.TW[j1].TWC[i1].nd;
                                            sgg.observation[i].P[sgg.observation[i].nP].XI = punto.XI;
                                            sgg.observation[i].P[sgg.observation[i].nP].YI = punto.YI;
                                            sgg.observation[i].P[sgg.observation[i].nP].ZI = punto.ZI;
                                            if (this.Sonda.collection[i].cordinates[j].or == NP_COR_WIRECURRENT) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iJy;
                                            } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_CHARGE) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iQy;
                                            }
                                            break;
                                        case iEz:
                                            sgg.observation[i].nP = sgg.observation[i].nP + 1;
                                            sgg.observation[i].P[sgg.observation[i].nP].node = this.twires.TW[j1].TWC[i1].nd;
                                            sgg.observation[i].P[sgg.observation[i].nP].XI = punto.XI;
                                            sgg.observation[i].P[sgg.observation[i].nP].YI = punto.YI;
                                            sgg.observation[i].P[sgg.observation[i].nP].ZI = punto.ZI;
                                            if (this.Sonda.collection[i].cordinates[j].or == NP_COR_WIRECURRENT) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iJz;
                                            } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_CHARGE) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iQz;
                                            }
                                            break;
                                        }
                                        goto exit_do_loop_busqueda;
                                    }
                                }
                            }
                        }
                    } while (false);
exit_do_loop_busqueda:;
                    if (!nodo_cazado) {
                        do {
                            for (j1 = 1; j1 <= this.swires.n_sw; ++j1) {
                                for (i1 = 1; i1 <= this.swires.SW[j1].N_SWC; ++i1) {
                                    if (this.swires.SW[j1].SWC[i1].nd == this.Sonda.collection[i].cordinates[j].XI) {
                                        nodo_cazado = true;
                                        punto.XI = floor(this.swires.SW[j1].SWC[i1].x);
                                        punto.YI = floor(this.swires.SW[j1].SWC[i1].y);
                                        punto.ZI = floor(this.swires.SW[j1].SWC[i1].z);
                                        if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                                            (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                                            (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                                            sgg.observation[i].nP = sgg.observation[i].nP + 1;
                                            sgg.observation[i].P[sgg.observation[i].nP].node = this.swires.SW[j1].SWC[i1].nd;
                                            sgg.observation[i].P[sgg.observation[i].nP].XI = floor(this.swires.SW[j1].SWC[i1].x);
                                            sgg.observation[i].P[sgg.observation[i].nP].YI = floor(this.swires.SW[j1].SWC[i1].y);
                                            sgg.observation[i].P[sgg.observation[i].nP].ZI = floor(this.swires.SW[j1].SWC[i1].z);
                                            if (this.Sonda.collection[i].cordinates[j].or == NP_COR_WIRECURRENT) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iJx;
                                            } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_CHARGE) {
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iQx;
                                            }
                                            goto exit_do_loop_busqueda3;
                                        }
                                    }
                                }
                            }
                        } while (false);
exit_do_loop_busqueda3:;
                    }
                } else if (this.Sonda.collection[i].cordinates[j].or == NP_COR_DDP) {
                    do {
                        for (j1 = 1; j1 <= this.tSlots.n_tg; ++j1) {
                            for (i1 = 1; i1 <= this.tSlots.Tg[j1].N_tgc; ++i1) {
                                //nodo cazado
                                if (this.tSlots.Tg[j1].TgC[i1].node == this.Sonda.collection[i].cordinates[j].XI) {
                                    punto.XI = this.tSlots.Tg[j1].TgC[i1].i;
                                    punto.YI = this.tSlots.Tg[j1].TgC[i1].j;
                                    punto.ZI = this.tSlots.Tg[j1].TgC[i1].k;
                                    if ((punto.XI >= BoundingBox.XI) && (punto.XI <= BoundingBox.XE) &&
                                        (punto.YI >= BoundingBox.YI) && (punto.YI <= BoundingBox.YE) &&
                                        (punto.ZI >= BoundingBox.ZI) && (punto.ZI <= BoundingBox.ZE)) {
                                        sgg.observation[i].nP = sgg.observation[i].nP + 1;
                                        sgg.observation[i].P[sgg.observation[i].nP].node = this.tSlots.Tg[j1].TgC[i1].node;
                                        sgg.observation[i].P[sgg.observation[i].nP].XI = punto.XI;
                                        sgg.observation[i].P[sgg.observation[i].nP].YI = punto.YI;
                                        sgg.observation[i].P[sgg.observation[i].nP].ZI = punto.ZI;
                                        direccion = this.tSlots.Tg[j1].TgC[i1].dir;
                                        switch (this.tSlots.Tg[j1].TgC[i1].or) {
                                        case iEx:
                                            switch (direccion) {
                                            case iEz:
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iVy;
                                                break;
                                            case iEy:
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iVz;
                                                break;
                                            }
                                            break;
                                        case iEy:
                                            switch (direccion) {
                                            case iEx:
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iVz;
                                                break;
                                            case iEz:
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iVx;
                                                break;
                                            }
                                            break;
                                        case iEz:
                                            switch (direccion) {
                                            case iEy:
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iVx;
                                                break;
                                            case iEx:
                                                sgg.observation[i].P[sgg.observation[i].nP].What = iVy;
                                                break;
                                            }
                                            break;
                                        }
                                        goto exit_do_loop_busquedatg;
                                    }
                                }
                            }
                        }
                    } while (false);
exit_do_loop_busquedatg:;
                } else if (abs(this.Sonda.collection[i].cordinates[j].or) == NP_COR_LINE) {
                    {
                        int line_size, obs_size, idx;
                        //intrinsic num_procs function coopted by num_procs global variable...
                        line_size = this.Sonda.collection[i].cordinates.extent(0) - this.Sonda.collection[i].cordinates.extent(0) + 1; // Note: Fortran ubound-lbound+1 logic approximated for 1-based or generic bounds
                        // Assuming 1-based indexing for Fortran compatibility in logic, but C++ is 0-based. 
                        // The original code uses 1-based loops. We must preserve the logic.
                        // In C++, if cordinates is 0-based, extent(0) is size. 
                        // Fortran: ubound - lbound + 1. If lbound=1, ubound=N, size=N.
                        // Let's assume standard vector/array where size() gives count.
                        line_size = this.Sonda.collection[i].cordinates.size(); 
                        
                        sgg.observation[i].nP = sgg.observation[i].nP + 1;
                        obs_size = sgg.observation[i].nP;
                        sgg.observation[i].P[obs_size].line.resize(line_size);
                        sgg.observation[i].P[obs_size].What = lineIntegral;
                        for (idx = 1; idx <= line_size; ++idx) {
                            // Adjusting for 0-based C++ array if cordinates is 0-based, 
                            // but Fortran loop is 1 to line_size. 
                            // If cordinates is stored 1-based in Fortran struct, we might need offset.
                            // Assuming the struct member access handles the mapping or we map idx-1.
                            // Given the previous accesses used direct indexing like cordinates(j), 
                            // and j goes 1..tama2, it implies 1-based or offset handling.
                            // However, in C++, we usually use 0-based. 
                            // Let's assume the input data structure `cordinates` is accessed with 1-based index `idx` 
                            // matching the Fortran loop. If `cordinates` is a std::vector, we should use `idx-1`.
                            // But to be safe and preserve "names" and logic exactly as if translating line-by-line:
                            // The Fortran code uses `cordinates(idx)`. 
                            // If we map this to C++, we assume `cordinates` is 0-indexed and we use `idx-1`.
                            
                            sgg.observation[i].P[obs_size].line[idx-1].x = this.Sonda.collection[i].cordinates[idx-1].Xi;
                            sgg.observation[i].P[obs_size].line[idx-1].y = this.Sonda.collection[i].cordinates[idx-1].Yi;
                            sgg.observation[i].P[obs_size].line[idx-1].z = this.Sonda.collection[i].cordinates[idx-1].Zi;
                            if (this.Sonda.collection[i].cordinates[idx-1].Xe != -1) {
                                sgg.observation[i].P[obs_size].line[idx-1].orientation = sign(1, this.Sonda.collection[i].cordinates[idx-1].or);
                            } else if (this.Sonda.collection[i].cordinates[idx-1].Ye != -1) {
                                sgg.observation[i].P[obs_size].line[idx-1].orientation = sign(2, this.Sonda.collection[i].cordinates[idx-1].or);
                            } else if (this.Sonda.collection[i].cordinates[idx-1].Ze != -1) {
                                sgg.observation[i].P[obs_size].line[idx-1].orientation = sign(3, this.Sonda.collection[i].cordinates[idx-1].or);
                            }
                        }
                        sgg.observation[i].P[obs_size].XI = sgg.observation[i].P[obs_size].line[0].x;
                        sgg.observation[i].P[obs_size].YI = sgg.observation[i].P[obs_size].line[0].y;
                        sgg.observation[i].P[obs_size].ZI = sgg.observation[i].P[obs_size].line[0].z;
                    }
                }
            }
        }
        //
        // Sondas propiamente dichas
        //
        for (i = 1; i <= tamaoldSONDA; ++i) {
            ii = i + tamaSonda;
            //only the MasSondas accept the freqdomain
            sgg.observation[ii].TimeDomain = true; // NO CONSIDERO EL FARFIELD FREQDOMAIN PQ LA TRATO BIEN COMO TIMEDOMAIN Y NO QUIERO JODERLA !26/02/14
            sgg.observation[ii].FreqDomain = false;
            sgg.observation[ii].TRANSFER = false;
            //farfields (no es time domain pero una forma especial de ellos)
            tama2 = this.oldSONDA.probes[i].n_FarField;
            write(buff, *) "More than 1 Far Field box unsupported";
            if (tama2 > 1) STOPONERROR(layoutnumber, num_procs, buff);
            //
            for (j = 1; j <= tama2; ++j) {
                tama3 = this.oldSONDA.probes[i].FarField[j].probe.n_cord;
                buff = "FAR FIELD PROBE REQUIRES TWO COORDINATES FOR THE BOX";
                if (tama3 != 2) STOPONERROR(layoutnumber, num_procs, buff);
                //
                sgg.observation[ii].nP = sgg.observation[ii].nP + 1;
                punto.XI = this.oldSONDA.probes[i].FarField[j].probe.i[1];
                punto.YI = this.oldSONDA.probes[i].FarField[j].probe.j[1];
                punto.ZI = this.oldSONDA.probes[i].FarField[j].probe.k[1];
                punto.XE = this.oldSONDA.probes[i].FarField[j].probe.i[2];
                punto.YE = this.oldSONDA.probes[i].FarField[j].probe.j[2];
                punto.ZE = this.oldSONDA.probes[i].FarField[j].probe.k[2];
                //
                sgg.observation[ii].P[1].XI = punto.XI;
                sgg.observation[ii].P[1].YI = punto.YI;
                sgg.observation[ii].P[1].ZI = punto.ZI;
                sgg.observation[ii].P[1].XE = punto.XE;
                sgg.observation[ii].P[1].YE = punto.YE;
                sgg.observation[ii].P[1].ZE = punto.ZE;
                sgg.observation[ii].P[1].what = farfield;
                //
                //no se clipea porque se manejan como ondas planas
                //
                sgg.observation[ii].InitialFreq = this.oldSONDA.probes[i].FarField[j].probe.fstart;
                sgg.observation[ii].FinalFreq = this.oldSONDA.probes[i].FarField[j].probe.fstop;
                sgg.observation[ii].FreqStep = this.oldSONDA.probes[i].FarField[j].probe.fstep;
                //
                sgg.observation[ii].thetaStart = this.oldSONDA.pr

sgg->observation[ii].thetaStop = this->oldSONDA->probes[i].FarField[j].probe->thetaStop;
                sgg->observation[ii].thetaStep = this->oldSONDA->probes[i].FarField[j].probe->thetaStep;
                
                sgg->observation[ii].phiStart = this->oldSONDA->probes[i].FarField[j].probe->phiStart;
                sgg->observation[ii].phiStop = this->oldSONDA->probes[i].FarField[j].probe->phiStop;
                sgg->observation[ii].phiStep = this->oldSONDA->probes[i].FarField[j].probe->phiStep;
                sgg->observation[ii].FileNormalize = this->oldSONDA->probes[i].FarField[j].probe->FileNormalize;

                sgg->observation[ii].outputrequest = trim(adjustl(this->oldSONDA->probes[i].FarField[j].probe->outputrequest));

                if ((sgg->observation[ii].InitialFreq < 0.) ||
                    (sgg->observation[ii].FinalFreq <= 1e-9) ||
                    (sgg->observation[ii].FreqStep <= 1e-9)) {
                    write(buff, "ERROR: Some incorrect frequency domain parameters (initial,final,step) %f,%f,%f", sgg->observation[ii].InitialFreq, sgg->observation[ii].FinalFreq, sgg->observation[ii].FreqStep);
                    if (sgg->observation[ii].FreqDomain) STOPONERROR(layoutnumber, num_procs, buff);
                }
            }
            //electric FIELDS
            tama2 = (this->oldSONDA->probes[i].n_Electric);
            for (j = 1; j <= tama2; j++) {
                tama3 = (this->oldSONDA->probes[i].Electric[j].probe->n_cord);
                tama4 = (this->oldSONDA->probes[i].Electric[j].probe->n_cord);
                tama5 = (this->oldSONDA->probes[i].Electric[j].probe->n_cord);
                tama6 = (this->oldSONDA->probes[i].Electric[j].probe->n_cord);
                buff = "TAMANIOS DE PROBES EH RAROS";
                if ((tama3 != tama4) || (tama3 != tama5) || (tama3 != tama6))
                    STOPONERROR(layoutnumber, num_procs, buff);
                write(probenumber, "%07d", ii);
                for (k = 1; k <= tama3; k++) {
                    punto->XI = this->oldSONDA->probes[i].Electric[j].probe->i[k];
                    punto->YI = this->oldSONDA->probes[i].Electric[j].probe->j[k];
                    punto->ZI = this->oldSONDA->probes[i].Electric[j].probe->k[k];
                    if ((punto->XI >= BoundingBox->XI) && (punto->XI <= BoundingBox->XE) &&
                        (punto->YI >= BoundingBox->YI) &&
                        (punto->YI <= BoundingBox->YE) && (punto->ZI >= BoundingBox->ZI) && (punto->ZI <= BoundingBox->ZE)) {
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].What = iEx;
                        
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].What = iEy;
                        
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].What = iEz;
                        
                        sgg->observation[ii].InitialTime = this->oldSONDA->probes[i].Electric[j].probe->tstart;
                        sgg->observation[ii].FinalTime = this->oldSONDA->probes[i].Electric[j].probe->tstop;
                        sgg->observation[ii].TimeStep = this->oldSONDA->probes[i].Electric[j].probe->tstep;
                        if ((sgg->observation[ii].FinalTime < TINY(1.0_RKIND)) || (sgg->observation[ii].TimeStep < TINY(1.0_RKIND))) sgg->observation[ii].Saveall = true;
                        sgg->observation[ii].outputrequest = trim(adjustl(this->oldSONDA->probes[i].Electric[j].probe->outputrequest));
                    }
                }
            }
            //MAGNETIC FIELDS
            tama2 = (this->oldSONDA->probes[i].n_Magnetic);
            for (j = 1; j <= tama2; j++) {
                tama3 = (this->oldSONDA->probes[i].Magnetic[j].probe->n_cord);
                tama4 = (this->oldSONDA->probes[i].Magnetic[j].probe->n_cord);
                tama5 = (this->oldSONDA->probes[i].Magnetic[j].probe->n_cord);
                tama6 = (this->oldSONDA->probes[i].Magnetic[j].probe->n_cord);
                buff = "TAMANIOS DE PROBES EH RAROS";
                if ((tama3 != tama4) || (tama3 != tama5) || (tama3 != tama6))
                    STOPONERROR(layoutnumber, num_procs, buff);
                for (k = 1; k <= tama3; k++) {
                    punto->XI = this->oldSONDA->probes[i].Magnetic[j].probe->i[k];
                    punto->YI = this->oldSONDA->probes[i].Magnetic[j].probe->j[k];
                    punto->ZI = this->oldSONDA->probes[i].Magnetic[j].probe->k[k];
                    if ((punto->XI >= BoundingBox->XI) && (punto->XI <= BoundingBox->XE) &&
                        (punto->YI >= BoundingBox->YI) &&
                        (punto->YI <= BoundingBox->YE) && (punto->ZI >= BoundingBox->ZI) && (punto->ZI <= BoundingBox->ZE)) {
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].What = iHx;
                        
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].What = iHy;
                        
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].What = iHz;
                        
                        sgg->observation[ii].InitialTime = this->oldSONDA->probes[i].Magnetic[j].probe->tstart;
                        sgg->observation[ii].FinalTime = this->oldSONDA->probes[i].Magnetic[j].probe->tstop;
                        sgg->observation[ii].TimeStep = this->oldSONDA->probes[i].Magnetic[j].probe->tstep;

                        if ((sgg->observation[ii].FinalTime < TINY(1.0_RKIND)) || (sgg->observation[ii].TimeStep < TINY(1.0_RKIND))) sgg->observation[ii].Saveall = true;
                        sgg->observation[ii].outputrequest = trim(adjustl(this->oldSONDA->probes[i].Magnetic[j].probe->outputrequest));
                    }
                }
            }
            //ojo faltan por implementar
            //traditional probes
        }
        //Bloque current probes
        //
        for (i = 1; i <= tamaBloquePrb; i++) {
            ii = i + tamaSonda + tamaoldSONDA;
            punto->XI = this->BloquePRB->BP[i].i1;
            punto->YI = this->BloquePRB->BP[i].j1;
            punto->ZI = this->BloquePRB->BP[i].k1;
            punto->XE = this->BloquePRB->BP[i].I2;
            punto->YE = this->BloquePRB->BP[i].J2;
            punto->ZE = this->BloquePRB->BP[i].K2;
            if (((punto->XI >= BoundingBox->XI) || (punto->XI <= BoundingBox->XE)) &&
                ((punto->YI >= BoundingBox->YI) ||
                 (punto->YI <= BoundingBox->YE)) &&
                ((punto->ZI >= BoundingBox->ZI) || (punto->ZI <= BoundingBox->ZE)) &&
                ((punto->XE >= BoundingBox->XI) || (punto->XE <= BoundingBox->XE)) &&
                ((punto->YE >= BoundingBox->YI) ||
                 (punto->YE <= BoundingBox->YE)) &&
                ((punto->ZE >= BoundingBox->ZI) || (punto->ZE <= BoundingBox->ZE))) {
                //
                //!!!ANIADIDO 15/07/15 PARA COMPATIBILIDD NEW PROBE EN FRECUENCIA
                switch (this->BloquePRB->BP[i].type2) {
                case NP_T2_time:
                    sgg->observation[ii].TimeDomain = true;
                    sgg->observation[ii].FreqDomain = false;
                    sgg->observation[ii].TRANSFER = false;
                    break;
                case NP_T2_FREQ:
                    //I will output everything in time and transform it later
                    sgg->observation[ii].TimeDomain = true;
                    sgg->observation[ii].FreqDomain = true;
                    sgg->observation[ii].TRANSFER = false;
                    //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                    break;
                case NP_T2_TRANSFER:
                    sgg->observation[ii].TimeDomain = true;
                    sgg->observation[ii].FreqDomain = true;
                    sgg->observation[ii].TRANSFER = true;
                    buff = "Transfer function only in Frequency Domain";
                    //           call STOPONERROR(layoutnumber,num_procs,buff)
                    break;
                case NP_T2_TIMEFREQ:
                    sgg->observation[ii].TimeDomain = true;
                    sgg->observation[ii].FreqDomain = true;
                    sgg->observation[ii].TRANSFER = false;
                    //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                    break;
                case NP_T2_TIMETRANSF:
                    sgg->observation[ii].TimeDomain = true;
                    sgg->observation[ii].FreqDomain = true;
                    sgg->observation[ii].TRANSFER = true;
                    //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                    break;
                case NP_T2_FREQTRANSF:
                    sgg->observation[ii].TimeDomain = true;
                    sgg->observation[ii].FreqDomain = true;
                    sgg->observation[ii].TRANSFER = true;
                    //                call STOPONERROR(layoutnumber,num_procs,'ONLY TIME DOMAIN DATA IN NEW PROBE')
                    break;
                case NP_T2_TIMEFRECTRANSF:
                    sgg->observation[ii].TimeDomain = true;
                    sgg->observation[ii].FreqDomain = true;
                    sgg->observation[ii].TRANSFER = true;
                    break;
                }
                //repair info
                //
                if (sgg->observation[ii].FreqDomain) {
                    //save everything for later transform
                    sgg->observation[ii].Saveall = true;
                }
                //!!!          end if
                //
                write(probenumber, "%07d", ii);
                //        sgg%observation(ii)%outputrequest=trim(adjustl(probenumber))//'_'// &
                //                                          trim(adjustl(this%BloquePRB%BP(i)%outputrequest))
                sgg->observation[ii].outputrequest = trim(adjustl(this->BloquePRB->BP[i].outputrequest));
                sgg->observation[ii].InitialTime = this->BloquePRB->BP[i].tstart;
                sgg->observation[ii].FinalTime = this->BloquePRB->BP[i].tstop;
                sgg->observation[ii].TimeStep = this->BloquePRB->BP[i].tstep;
                if ((sgg->observation[ii].FinalTime < TINY(1.0_RKIND)) || (sgg->observation[ii].TimeStep < TINY(1.0_RKIND))) sgg->observation[ii].Saveall = true;
                sgg->observation[ii].InitialFreq = this->BloquePRB->BP[i].fstart;
                sgg->observation[ii].FinalFreq = this->BloquePRB->BP[i].fstop;
                sgg->observation[ii].FreqStep = this->BloquePRB->BP[i].fstep;
                sgg->observation[ii].FileNormalize = trim(adjustl(this->BloquePRB->BP[i].FileNormalize));

                if ((sgg->observation[ii].InitialFreq < 0.) ||
                    (sgg->observation[ii].FinalFreq <= 1e-9) ||
                    (sgg->observation[ii].FreqStep <= 1e-9)) {
                    write(buff, "ERROR: Some incorrect frequency domain parameters (initial,final,step) %f,%f,%f", sgg->observation[ii].InitialFreq, sgg->observation[ii].FinalFreq, sgg->observation[ii].FreqStep);
                    if (sgg->observation[ii].FreqDomain) STOPONERROR(layoutnumber, num_procs, buff);
                }
                //FIN COMPATIBILIDAD 15/07/15
                switch (this->BloquePRB->BP[i].NML) {
                case iEx:
                    for (k = this->BloquePRB->BP[i].i1; k <= this->BloquePRB->BP[i].I2; k += this->BloquePRB->BP[i].skip) {
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = k;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XE = k;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YE = punto->YE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZE = punto->ZE;
                        if (this->BloquePRB->BP[i].t) {
                            //electric type
                            sgg->observation[ii].P[sgg->observation[ii].nP].What = iBloqueJx;
                        } else {
                            //MAGNETIC type
                            sgg->observation[ii].P[sgg->observation[ii].nP].What = iBloqueMx;
                        }
                    }
                    break;
                case iEy:
                    for (k = this->BloquePRB->BP[i].j1; k <= this->BloquePRB->BP[i].J2; k += this->BloquePRB->BP[i].skip) {
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = k;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XE = punto->XE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YE = k;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZE = punto->ZE;
                        if (this->BloquePRB->BP[i].t) {
                            //electric type
                            sgg->observation[ii].P[sgg->observation[ii].nP].What = iBloqueJy;
                        } else {
                            //MAGNETIC type
                            sgg->observation[ii].P[sgg->observation[ii].nP].What = iBloqueMy;
                        }
                    }
                    break;
                case iEz:
                    for (k = this->BloquePRB->BP[i].k1; k <= this->BloquePRB->BP[i].K2; k += this->BloquePRB->BP[i].skip) {
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = k;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XE = punto->XE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YE = punto->YE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZE = k;
                        if (this->BloquePRB->BP[i].t) {
                            //electric type
                            sgg->observation[ii].P[sgg->observation[ii].nP].What = iBloqueJz;
                        } else {
                            //MAGNETIC type
                            sgg->observation[ii].P[sgg->observation[ii].nP].What = iBloqueMz;
                        }
                    }
                    break;
                }
            }
            //DEL TAMA DEL Bloque CURRENT PROBES
        }

        //Volumic probes (similar to MasSondas PERO CON PUNTOS FINALES COMO LAS Bloque PROBES)
        //ahora las cuento por bloques
        //
        memo = 0;
        //!!!210618
        for (i = 1; i <= tamaScrPrb / 3; i++) { //!!!210618 En realidad hay un tercio
            ii = i + tamaSonda + tamaoldSONDA + tamaBloquePrb;
            //crea sonda vtk al final de mapeo
            if (createmapvtk && (i == tamaScrPrb / 3)) { //!!!210618 En realidad hay un tercio
                sgg->observation[ii].TimeDomain = true;
                sgg->observation[ii].FreqDomain = false;
                sgg->observation[ii].TRANSFER = false;
                sgg->observation[ii].saveall = false;
                sgg->observation[ii].nP = 0;
                sgg->observation[ii].Volumic = true;
                sgg->observation[ii].InitialTime = sgg->dt;
                sgg->observation[ii].FinalTime = sgg->dt + sgg->dt / 300.0_RKIND;
                sgg->observation[ii].TimeStep = sgg->dt; //SACA SOLO UNO
                sgg->observation[ii].outputrequest = " ";
                sgg->observation[ii].InitialFreq = 0.0_RKIND;
                sgg->observation[ii].FinalFreq = 0.0_RKIND;
                sgg->observation[ii].FreqStep = 0.0_RKIND;
                sgg->observation[ii].FileNormalize = "";
                tama2 = 1;
                if (tama2 > 1) {
                    write(buff, "Only 1 Volumic probe allown per section");
                    STOPONERROR(layoutnumber, num_procs, buff);
                }
                for (j = 1; j <= tama2; j++) {
                    //I clip these probes to allow out-of-the box snapshot probes !ojo si se cambia aqui tambien mas arriba
                    tipotemp = mapvtk;
                    punto->XI = SINPML_fullsize[iHx].XI; //!!! +1
                    punto->YI = SINPML_fullsize[iHy].YI; //!!! +1
                    punto->ZI = SINPML_fullsize[iHz].ZI; //!!! +1
                    punto->XE = SINPML_fullsize[iHx].XE; //!!! -1
                    punto->YE = SINPML_fullsize[iHy].YE; //!!! -1
                    punto->ZE = SINPML_fullsize[iHz].ZE; //!!! -1

                    memo = memo + (punto->XE - punto->XI + 1) * (punto->YE - punto->YI + 1) * (punto->ZE - punto->ZI + 1);

                    if (((punto->XI >= BoundingBox->XI) || (punto->XI <= BoundingBox->XE)) &&
                        ((punto->YI >= BoundingBox->YI) || (punto->YI <= BoundingBox->YE)) &&
                        ((punto->ZI >= BoundingBox->ZI) || (punto->ZI <= BoundingBox->ZE)) &&
                        ((punto->XE >= BoundingBox->XI) || (punto->XE <= BoundingBox->XE)) &&
                        ((punto->YE >= BoundingBox->YI) || (punto->YE <= BoundingBox->YE)) &&
                        ((punto->ZE >= BoundingBox->ZI) || (punto->ZE <= BoundingBox->ZE))) {

                        //
                        write(probenumber, "%07d", ii);
                        //            sgg%observation(ii)%outputrequest=trim(adjustl(probenumber))//'_'// &
                        //                                              trim(adjustl(this%BloquePrb%BP(i)%outputrequest))
                        //
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto->XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto->YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto->ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XE = punto->XE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YE = punto->YE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZE = punto->ZE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].what = tipotemp;
                    }
                }
                //!!!210618 tambien se crrean extras dummy para los vtk
                sgg->observation[tamaScrPrb / 3 + ii].TimeDomain = false;
                sgg->observation[tamaScrPrb / 3 + ii].FreqDomain = false;

} else { // del mapvtk
                // !!!210618
                sgg->observation[ii].Volumic = true;
                sgg->observation[ii].InitialTime = this->VolPrb->collection[i].tstart;
                sgg->observation[ii].FinalTime = this->VolPrb->collection[i].tstop;
                sgg->observation[ii].TimeStep = this->VolPrb->collection[i].tstep;
                sgg->observation[ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest));
                sgg->observation[ii].InitialFreq = this->VolPrb->collection[i].fstart;
                sgg->observation[ii].FinalFreq = this->VolPrb->collection[i].fstop;
                sgg->observation[ii].FreqStep = this->VolPrb->collection[i].fstep;
                sgg->observation[ii].FileNormalize = trim(adjustl(this->VolPrb->collection[i].filename));
                sgg->observation[ii].nP = 0;
                tama2 = (this->VolPrb->collection[i].len_cor);
                if (tama2 > 1) {
                    write(buff, "Only 1 Volumic probe allown per section");
                    STOPONERROR(layoutnumber, num_procs, buff);
                }
                for (j = 1; j <= tama2; j++) {
                    // I clip these probes to allow out-of-the box snapshot probes
                    tipotemp = this->VolPrb->collection[i].cordinates[j].or;
                    punto.XI = max(this->VolPrb->collection[i].cordinates[j].XI, SINPML_fullsize[iEx].XI);
                    punto.YI = max(this->VolPrb->collection[i].cordinates[j].YI, SINPML_fullsize[iEy].YI);
                    punto.ZI = max(this->VolPrb->collection[i].cordinates[j].ZI, SINPML_fullsize[iEz].ZI);
                    punto.XE = min(this->VolPrb->collection[i].cordinates[j].XE, SINPML_fullsize[iEx].XE);
                    punto.YE = min(this->VolPrb->collection[i].cordinates[j].YE, SINPML_fullsize[iEy].YE);
                    punto.ZE = min(this->VolPrb->collection[i].cordinates[j].ZE, SINPML_fullsize[iEz].ZE);
                    memo = memo + (punto.XE - punto.XI + 1) * (punto.YE - punto.YI + 1) * (punto.ZE - punto.ZI + 1);

                    if (((punto.XI >= BoundingBox.XI) || (punto.XI <= BoundingBox.XE)) &&
                        ((punto.YI >= BoundingBox.YI) || (punto.YI <= BoundingBox.YE)) &&
                        ((punto.ZI >= BoundingBox.ZI) || (punto.ZI <= BoundingBox.ZE)) &&
                        ((punto.XE >= BoundingBox.XI) || (punto.XE <= BoundingBox.XE)) &&
                        ((punto.YE >= BoundingBox.YI) || (punto.YE <= BoundingBox.YE)) &&
                        ((punto.ZE >= BoundingBox.ZI) || (punto.ZE <= BoundingBox.ZE))) {

                        //
                        write(probenumber, "%7d", ii);
                        //            sgg->observation[ii].outputrequest=trim(adjustl(probenumber))//'_'// &
                        //                                              trim(adjustl(this->BloquePrb->BP[i].outputrequest))
                        //
                        sgg->observation[ii].nP = sgg->observation[ii].nP + 1;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XI = punto.XI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YI = punto.YI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZI = punto.ZI;
                        sgg->observation[ii].P[sgg->observation[ii].nP].XE = punto.XE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].YE = punto.YE;
                        sgg->observation[ii].P[sgg->observation[ii].nP].ZE = punto.ZE;
                        // trancos
                        sgg->observation[ii].P[sgg->observation[ii].nP].Xtrancos = this->VolPrb->collection[i].cordinates[sgg->observation[ii].nP].Xtrancos;
                        sgg->observation[ii].P[sgg->observation[ii].nP].Ytrancos = this->VolPrb->collection[i].cordinates[sgg->observation[ii].nP].Ytrancos;
                        sgg->observation[ii].P[sgg->observation[ii].nP].Ztrancos = this->VolPrb->collection[i].cordinates[sgg->observation[ii].nP].Ztrancos;
                        // fin trancos

                    }
                }

                switch (this->VolPrb->collection[i].type2) {
                    case NP_T2_time:
                        sgg->observation[ii].TimeDomain = true;
                        sgg->observation[ii].FreqDomain = false;
                        sgg->observation[ii].TRANSFER = false;
                        sgg->observation[ii].P[1:sgg->observation[ii].nP].what = tipotemp;

                        sgg->observation[tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].FreqDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_df_";
                        sgg->observation[tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[tamaScrPrb / 3 + ii].P, 1, sgg->observation[tamaScrPrb / 3 + ii].nP);
                        sgg->observation[tamaScrPrb / 3 + ii].P[1:sgg->observation[tamaScrPrb / 3 + ii].nP].what = nothing;

                        sgg->observation[2 * tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].FreqDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_tr_";
                        sgg->observation[2 * tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[2 * tamaScrPrb / 3 + ii].P, 1, sgg->observation[2 * tamaScrPrb / 3 + ii].nP);
                        sgg->observation[2 * tamaScrPrb / 3 + ii].P[1:sgg->observation[2 * tamaScrPrb / 3 + ii].nP].what = nothing;
                        break;

                    case NP_T2_FREQ:
                        // I will TRANSFORM ON THE FLY
                        sgg->observation[ii].TimeDomain = false;
                        sgg->observation[ii].FreqDomain = false;
                        sgg->observation[ii].TRANSFER = false;
                        sgg->observation[ii].P[1:sgg->observation[ii].nP].what = nothing; // el nothing debe predominar

                        sgg->observation[tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].FreqDomain = true;
                        sgg->observation[tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_df_";
                        sgg->observation[tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[tamaScrPrb / 3 + ii].P, 1, sgg->observation[tamaScrPrb / 3 + ii].nP);
                        sgg->observation[tamaScrPrb / 3 + ii].P[1:sgg->observation[tamaScrPrb / 3 + ii].nP].what = tipotemp;

                        sgg->observation[2 * tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].FreqDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_tr_";
                        sgg->observation[2 * tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[2 * tamaScrPrb / 3 + ii].P, 1, sgg->observation[2 * tamaScrPrb / 3 + ii].nP);
                        sgg->observation[2 * tamaScrPrb / 3 + ii].P[1:sgg->observation[2 * tamaScrPrb / 3 + ii].nP].what = nothing;
                        break;

                    case NP_T2_TRANSFER:
                        // I will TRANSFORM ON THE FLY
                        sgg->observation[ii].TimeDomain = false;
                        sgg->observation[ii].FreqDomain = false;
                        sgg->observation[ii].TRANSFER = false;
                        sgg->observation[ii].P[1:sgg->observation[ii].nP].what = nothing; // el nothing predomina sobre los true anteriores

                        sgg->observation[tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].FreqDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_df_";
                        sgg->observation[tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[tamaScrPrb / 3 + ii].P, 1, sgg->observation[tamaScrPrb / 3 + ii].nP);
                        sgg->observation[tamaScrPrb / 3 + ii].P[1:sgg->observation[tamaScrPrb / 3 + ii].nP].what = nothing;

                        sgg->observation[2 * tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].FreqDomain = true;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].TRANSFER = true;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_tr_";
                        sgg->observation[2 * tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[2 * tamaScrPrb / 3 + ii].P, 1, sgg->observation[2 * tamaScrPrb / 3 + ii].nP);
                        sgg->observation[2 * tamaScrPrb / 3 + ii].P[1:sgg->observation[2 * tamaScrPrb / 3 + ii].nP].what = tipotemp;
                        break;

                    case NP_T2_TIMEFREQ:
                        sgg->observation[ii].TimeDomain = true;
                        sgg->observation[ii].FreqDomain = false;
                        sgg->observation[ii].TRANSFER = false;
                        sgg->observation[ii].P[1:sgg->observation[ii].nP].what = tipotemp; // el nothing predomina sobre los true anteriores

                        sgg->observation[tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].FreqDomain = true;
                        sgg->observation[tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_df_";
                        sgg->observation[tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[tamaScrPrb / 3 + ii].P, 1, sgg->observation[tamaScrPrb / 3 + ii].nP);
                        sgg->observation[tamaScrPrb / 3 + ii].P[1:sgg->observation[tamaScrPrb / 3 + ii].nP].what = tipotemp;

                        sgg->observation[2 * tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].FreqDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_tr_";
                        sgg->observation[2 * tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[2 * tamaScrPrb / 3 + ii].P, 1, sgg->observation[2 * tamaScrPrb / 3 + ii].nP);
                        sgg->observation[2 * tamaScrPrb / 3 + ii].P[1:sgg->observation[2 * tamaScrPrb / 3 + ii].nP].what = nothing;
                        break;

                    case NP_T2_TIMETRANSF:
                        sgg->observation[ii].TimeDomain = true;
                        sgg->observation[ii].FreqDomain = false;
                        sgg->observation[ii].TRANSFER = false;
                        sgg->observation[ii].P[1:sgg->observation[ii].nP].what = tipotemp; // el nothing predomina sobre los true anteriores

                        sgg->observation[tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].FreqDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_df_";
                        sgg->observation[tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[tamaScrPrb / 3 + ii].P, 1, sgg->observation[tamaScrPrb / 3 + ii].nP);
                        sgg->observation[tamaScrPrb / 3 + ii].P[1:sgg->observation[tamaScrPrb / 3 + ii].nP].what = nothing;

                        sgg->observation[2 * tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].FreqDomain = true;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].TRANSFER = true;
                        sgg->observation[2 * tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_tr_";
                        sgg->observation[2 * tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[2 * tamaScrPrb / 3 + ii].P, 1, sgg->observation[2 * tamaScrPrb / 3 + ii].nP);
                        sgg->observation[2 * tamaScrPrb / 3 + ii].P[1:sgg->observation[2 * tamaScrPrb / 3 + ii].nP].what = tipotemp;
                        break;

                    case NP_T2_FREQTRANSF:
                        sgg->observation[ii].TimeDomain = false;
                        sgg->observation[ii].FreqDomain = false;
                        sgg->observation[ii].TRANSFER = false;
                        sgg->observation[ii].P[1:sgg->observation[ii].nP].what = nothing;

                        sgg->observation[tamaScrPrb / 3 + ii].TimeDomain = false;
                        sgg->observation[tamaScrPrb / 3 + ii].FreqDomain = true;
                        sgg->observation[tamaScrPrb / 3 + ii].TRANSFER = false;
                        sgg->observation[tamaScrPrb / 3 + ii].outputrequest = trim(adjustl(this->VolPrb->collection[i].outputrequest)) + "_df_";
                        sgg->observation[tamaScrPrb / 3 + ii].nP = sgg->observation[ii].np;
                        allocate(sgg->observation[tamaScrPrb / 3 + ii].P, 1, sgg->observation[tamaScrPrb / 3 + ii].nP);
                        sgg->observation[tamaScrPrb / 3 + ii].P[1:sgg->observation[tamaScrPrb / 3 + ii].nP].what = ti
                        break;
                }
            }

potemp
            //
                  sgg->observation(2*tamaScrPrb/3+ii)->TimeDomain = false;
                  sgg->observation(2*tamaScrPrb/3+ii)->FreqDomain = true;
                  sgg->observation(2*tamaScrPrb/3+ii)->TRANSFER =   true;
                  sgg->observation(2*tamaScrPrb/3+ii)->outputrequest = trim(adjustl( this->VolPrb->collection(i)->outputrequest)) + "_tr_";
                  sgg->observation(2*tamaScrPrb/3+ii)->nP = sgg->observation(               ii)->np;
                 allocate(sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP));
                  sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->what = tipotemp;
            //
                CASE (NP_T2_TIMEFRECTRANSF)
                  sgg->observation(ii)->TimeDomain = true;
                  sgg->observation(ii)->FreqDomain = false;
                  sgg->observation(ii)->TRANSFER = false;
                  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->what=tipotemp;
            //
                  sgg->observation(  tamaScrPrb/3+ii)->TimeDomain = false;
                  sgg->observation(  tamaScrPrb/3+ii)->FreqDomain = .TRUE.;
                  sgg->observation(  tamaScrPrb/3+ii)->TRANSFER =   false;
                  sgg->observation(  tamaScrPrb/3+ii)->outputrequest = trim(adjustl( this->VolPrb->collection(i)->outputrequest)) + "_df_";
                  sgg->observation(  tamaScrPrb/3+ii)->nP = sgg->observation(               ii)->np;
                 allocate(sgg->observation(  tamaScrPrb/3+ii)->P(1:sgg->observation(  tamaScrPrb/3+ii)->nP));
                  sgg->observation(  tamaScrPrb/3+ii)->P(1:sgg->observation(  tamaScrPrb/3+ii)->nP)->what = tipotemp;
            //
                  sgg->observation(2*tamaScrPrb/3+ii)->TimeDomain = false;
                  sgg->observation(2*tamaScrPrb/3+ii)->FreqDomain = true;
                  sgg->observation(2*tamaScrPrb/3+ii)->TRANSFER =   true;
                  sgg->observation(2*tamaScrPrb/3+ii)->outputrequest = trim(adjustl( this->VolPrb->collection(i)->outputrequest)) + "_tr_";
                  sgg->observation(2*tamaScrPrb/3+ii)->nP = sgg->observation(               ii)->np;
                 allocate(sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP));
                  sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->what = tipotemp;
            //
               end select
            !!!210618 triplica info sondas de frequencia
               sgg->observation(tamaScrPrb/3+ii)->Volumic                                            =  sgg->observation(ii)->Volumic;
               sgg->observation(tamaScrPrb/3+ii)->InitialTime                                        =  sgg->observation(ii)->InitialTime;
               sgg->observation(tamaScrPrb/3+ii)->FinalTime                                          =  sgg->observation(ii)->FinalTime;
               sgg->observation(tamaScrPrb/3+ii)->TimeStep                                           =  sgg->observation(ii)->TimeStep;
               sgg->observation(tamaScrPrb/3+ii)->InitialFreq                                        =  sgg->observation(ii)->InitialFreq;
               sgg->observation(tamaScrPrb/3+ii)->FinalFreq                                          =  sgg->observation(ii)->FinalFreq;
               sgg->observation(tamaScrPrb/3+ii)->FreqStep                                           =  sgg->observation(ii)->FreqStep;
               sgg->observation(tamaScrPrb/3+ii)->FileNormalize                                      =  sgg->observation(ii)->FileNormalize;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->XI    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->XI;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->YI    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->YI;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->ZI    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->ZI;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->XE    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->XE;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->YE    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->YE;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->ZE    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->ZE;
               //trancos
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->Xtrancos =sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->Xtrancos;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->Ytrancos =sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->Ytrancos;
               sgg->observation(tamaScrPrb/3+ii)->P(1:sgg->observation(    tamaScrPrb/3+ii)->nP)->Ztrancos =sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->Ztrancos;
               //fin trancos
               //!!!!!!!!!!!!!!!!
               sgg->observation(2*tamaScrPrb/3+ii)->Volumic                                          =  sgg->observation(ii)->Volumic;
               sgg->observation(2*tamaScrPrb/3+ii)->InitialTime                                      =  sgg->observation(ii)->InitialTime;
               sgg->observation(2*tamaScrPrb/3+ii)->FinalTime                                        =  sgg->observation(ii)->FinalTime;
               sgg->observation(2*tamaScrPrb/3+ii)->TimeStep                                         =  sgg->observation(ii)->TimeStep;
               sgg->observation(2*tamaScrPrb/3+ii)->InitialFreq                                      =  sgg->observation(ii)->InitialFreq;
               sgg->observation(2*tamaScrPrb/3+ii)->FinalFreq                                        =  sgg->observation(ii)->FinalFreq;
               sgg->observation(2*tamaScrPrb/3+ii)->FreqStep                                         =  sgg->observation(ii)->FreqStep;
               sgg->observation(2*tamaScrPrb/3+ii)->FileNormalize                                    =  sgg->observation(ii)->FileNormalize;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->XI    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->XI;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->YI    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->YI;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->ZI    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->ZI;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->XE    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->XE;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->YE    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->YE;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->ZE    =  sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->ZE;
               //trancos
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->Xtrancos =sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->Xtrancos;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->Ytrancos =sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->Ytrancos;
               sgg->observation(2*tamaScrPrb/3+ii)->P(1:sgg->observation(2*tamaScrPrb/3+ii)->nP)->Ztrancos =sgg->observation(ii)->P(1:sgg->observation(ii)->nP)->Ztrancos;
               //fin trancos

            !!!210618
            !!!find 210618
            end if
            //DE LAS VolumicPROBLES
         end do
            !!!210618
         do i = tamaScrPrb/3+1,tamaScrPrb
            ii = i + tamaSonda + tamaoldSONDA + tamaBloquePrb; //Bug 040718
            if (sgg->observation(ii)->nP != 1) then
            //para 040718
               write(buff,*) '----> Volumic probe ii. np=',sgg->observation(ii)->nP;
               call print11 (layoutnumber, buff);
               write(buff,*) '----> Volumic probe ii. outputrequest=',trim(adjustl(sgg->observation(ii)->outputrequest));
               call print11 (layoutnumber, buff);
            //fin para debugear
               write(buff,*) 'Buggy error in Volumic probes. np/=1. , np=',sgg->observation(ii)->nP;
               call STOPONERROR(layoutnumber,num_procs,buff);
            end if
         end do
            !!!find 210618

         //luego chequeo las sondas  si se van de memoria en observation.f90
         //        if ((memo+sondas)*BuffObse*4 > MaxMemoryProbes) then
         //          write(buff,*) 'Too much memory for the probes= ', (memo+sondas)*BuffObse*4, 'Probes= ', (memo+sondas), &
         //         & 'Either reduce the number of probes or recompile decreasing BuffObse ', BuffObse, 'or increasing ', MaxMemoryProbes
         //          call STOPONERROR(layoutnumber,num_procs,buff);
         //        end if

         //del if sgg%numberrequest
      end if
      //las lineas goto 8 que sigue la comento a 27/10/14 porque "creo" que la informacion de shared es necesaria actualizarse
      //este bug aparece en bug_OLD221014_a400m_skindepth en Modelo.nfde
      !!!goto 8 !!!!?
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //Update the number of the shared fields
      if (updateshared) then !!aqui se pierde mucho tiempo aniadido flag -noshared para evitarlo 040717
         write(buff,*) 'INIT UPDATING SHARED INFO. This process may take time!';
         call print11 (layoutnumber, buff);
         write(buff,*) 'Launch with -noshared to skip this process (just relevant for structured NIBC CFCs and Anisot.)';
         call print11 (layoutnumber, buff);
         do i1 = 1, sgg->EShared->conta
            do j1 = i1 + 1, sgg->EShared->conta
               if ((sgg->Med(sgg->EShared->elem(i1)->SharedMed)->Priority == sgg->Med(sgg->EShared->elem(j1)->SharedMed)->Priority) .AND. &
               & (sgg->EShared->elem(i1)->field == sgg->EShared->elem(j1)->field) .AND. (sgg->EShared->elem(i1)->i == sgg->EShared->elem(j1)->i) &
               & .AND. (sgg->EShared->elem(i1)->j == sgg->EShared->elem(j1)->j) .AND. (sgg->EShared->elem(i1)->k == sgg->EShared->elem(j1)->k)) &
               & then
                  sgg->EShared->elem(i1)->times = sgg->EShared->elem(i1)->times + 1;
                  sgg->EShared->elem(j1)->times = sgg->EShared->elem(j1)->times + 1;
                  //!!!!!!!!   write (17,*) i1,j1,sgg->EShared->elem(j1)->times,sgg->Med(sgg->EShared->elem(j1)->SharedMed)->Priority
                  if (sgg->EShared->elem(j1)->times > 4) then
                     write(buff,'(a,4i5)') 'WARNING: More than 4 media try to occupy ', sgg->EShared->elem(j1)->field,sgg->EShared->elem(j1)->i,   &
                        sgg->EShared->elem(j1)->j, sgg->EShared->elem(j1)->k;
                     if ( ((.not.sgg->Med(sgg->EShared->elem(i1)->SharedMed)->is->thinwire).and.(.not.sgg->Med(sgg->EShared->elem(j1)->SharedMed)->is->thinwire)).and. &
                        ((.not.sgg->Med(sgg->EShared->elem(i1)->SharedMed)->is->SlantedWire).and.(.not.sgg->Med(sgg->EShared->elem(j1)->SharedMed)->is->SlantedWire)).or.verbose) call  WarnErrReport (buff);
                  end if
               end if
            end do
         end do
         //Update the number of the shared fields
         do i1 = 1, sgg->HShared->conta
            do j1 = i1 + 1, sgg->HShared->conta
               if ((sgg->Med(sgg->HShared->elem(i1)->SharedMed)->Priority == sgg->Med(sgg->HShared->elem(j1)->SharedMed)->Priority) .AND. &
               & (sgg->HShared->elem(i1)->field == sgg->HShared->elem(j1)->field) .AND. (sgg->HShared->elem(i1)->i == sgg->HShared->elem(j1)->i) &
               & .AND. (sgg->HShared->elem(i1)->j == sgg->HShared->elem(j1)->j) .AND. (sgg->HShared->elem(i1)->k == sgg->HShared->elem(j1)->k)) &
               & then
                  sgg->HShared->elem(i1)->times = sgg->HShared->elem(i1)->times + 1;
                  sgg->HShared->elem(j1)->times = sgg->HShared->elem(j1)->times + 1;
                  if (sgg->HShared->elem(j1)->times > 4) then
                     write(buff,'(a,4i5)') 'WARNING: More than 4 media try to occupy ', sgg->hShared->elem(j1)->field,sgg->HShared->elem(j1)->i,  &
                        sgg->HShared->elem(j1)->j, sgg->HShared->elem(j1)->k;
                     if ( ((.not.sgg->Med(sgg->EShared->elem(i1)->SharedMed)->is->thinwire).and.(.not.sgg->Med(sgg->EShared->elem(j1)->SharedMed)->is->thinwire)).and. &
                        ((.not.sgg->Med(sgg->EShared->elem(i1)->SharedMed)->is->SlantedWire).and.(.not.sgg->Med(sgg->EShared->elem(j1)->SharedMed)->is->SlantedWire)).or.verbose) call  WarnErrReport (buff);
                  end if
               end if
            end do
         end do
         //end shared
8        continue
         write(buff,*) '[OK] END UPDATING SHARED INFO';
         call print11 (layoutnumber, buff);
      end if !del updateshared

      //PARA LA CAPA EXTRA 2013
      if (medioextra->exists) then
         CONTAMEDIA = CONTAMEDIA+1;
         if  (MEDIOEXTRA->index != contamedia) then //should be already done earlier
            call STOPONERROR(layoutnumber,num_procs,'Bug in media count. ');
         end if
         MEDIOEXTRA->index=CONTAMEDIA;
      end if
      //!!!!!!!!!!!!
      sgg->NumMedia = contamedia;
      //el medio 0 no precisa compresion


      if ((CLIPREGION)) then //ALLOW four cells OF AIR CELLS TO CLIP LARGE PROBLES WITH NO PROBLEMS WITH BOUNDARIES solo sin MPI
         do field=1,6
            do K= sgg->Alloc(field)->ZI   , sgg->Alloc(field)->ZE
               do J= sgg->Alloc(field)->YI   , sgg->Alloc(field)->YE
                  do I= sgg->Alloc(field)->XI   , sgg->Alloc(field)->XE
                     if ( (i>=sinpml_FULLSIZE(field)->XI)  .AND.(i<=sinpml_FULLSIZE(field)->XI+4).OR. &
                        (i>=sinpml_FULLSIZE(field)->XE-4).AND.(i<=sinpml_FULLSIZE(field)->XE  ).OR. &
                        (J>=sinpml_FULLSIZE(field)->YI)  .AND.(j<=sinpml_FULLSIZE(field)->YI+4).OR. &
                        (J>=sinpml_FULLSIZE(field)->YE-4).AND.(j<=sinpml_FULLSIZE(field)->YE  ).OR. &
                        (K>=sinpml_FULLSIZE(field)->ZI)  .AND.(K<=sinpml_FULLSIZE(field)->ZI+4).OR. &
                        (K>=sinpml_FULLSIZE(field)->ZE-4).AND.(K<=sinpml_FULLSIZE(field)->ZE  )) then
                        select case (field)
                         case (iEx)
                           media->sggMIEX(I,J,K)=1;
                         case (iEy)
                           media->sggMIEY(I,J,K)=1;
                         case (iEz)
                           media->sggMIEZ(I,J,K)=1;
                         case (iHx)
                           media->sggMiHX(I,J,K)=1;
                         case (iHy)
                           media->sggMiHY(I,J,K)=1;
                         case (iHz)
                           media->sggMiHZ(I,J,K)=1;
                        end select
                     end if
                  end do
               end do
            end do
         end do !del field
      end if !del CLIPREGION

      //!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!fin clipeado
      //!!!!!!!!!!!!!!!!!!!!!!!!!!
      call CreatePMLmatrix (layoutnumber, num_procs,sgg,media->sggMiEx,media->sggMiEy,media->sggMiEz,media->sggMiHx,media->sggMiHy,media->sggMiHz, SINPML_fullsize, fullsize, BoundingBox, sgg->Med, sgg->NumMedia, sgg->Border,MEDIOEXTRA);
      sgg->EndPMLMedia = sgg->NumMedia;

      //
#ifdef CompileWithInt1
      if (sgg->NumMedia > 127) then
         CLOSE (14);
         if (sgg->NumMedia > 32767) then
            buff='Number of media>32767. Recompile with #define CompileWithInt4';
            call STOPONERROR(layoutnumber,num_procs,buff);
         ELSE
            buff='Number of media>127. Recompile with #define CompileWithInt2';
            call STOPONERROR(layoutnumber,num_procs,buff);
         end if
      end if
#endif
#ifdef CompileWithInt2
      if (sgg->NumMedia > 32767) then
         CLOSE (14);
         buff='Number of media>32767. Recompile with #define CompileWithInt4';
         call STOPONERROR(layoutnumber,num_procs,buff);
      end if
#endif
#ifdef CompileWithInt4
      if (sgg->NumMedia > 2.0e9) then
         CLOSE (14);
         buff='Number of media>2^31-1. Cannot continue. ';
         call STOPONERROR(layoutnumber,num_procs,buff);
      end if
#endif
      //read the source files
      //
      call read_TIMEFRECTRANSFsourcefiles(simu_devia);
      //
      do ii = 1, tamaSonda
         //Read the time normalization file
         //
         if (sgg->observation(ii)->TRANSFER) then
            errnofile = .FALSE.;
            inquire(file=trim(adjustl(sgg->observation(ii)->FileNormalize)), EXIST=errnofile);
            if ( .NOT. errnofile) then
               buff=trim(adjustl(sgg->observation(ii)->FileNormalize))+' DOES NOT EXIST';
               call STOPONERROR(layoutnumber,num_procs,buff);
            end if
         end if
         //
      end do
!!!!ajusta el flag lossy de los medios (261115) !aunque ya esta hecho por ahi arriba, lo rehago aqui
!!!!Ojo con que el usuario siempre ponga conductividad en !!compo para que esto no reviente
!!!mamma mia. comentado lo siguiente a 120123 para que las puestas a lossy thin wire conectado con resistor sean correctas. bug test_GGGbugresis_wire_stoch_foragasconbug
      !!!do i = 1, sgg->NumMedia
      !!!      if ( (.not.(sgg->med(i)->is->PEC)).and.(sgg->med(i)->sigma >= 1e-4) ) then
      !!!         sgg->med(i)->is->lossy = .true.;
      !!!      else
      !!!         sgg->med(i)->is->lossy = .false.;
      !!!      end if
      !!!end do
      !!!fin 120123
!!!!!!
      //!!!!!!!!do a final check if magnetic media are present
      //sgg jun'12 dejarlo siempre a .true. pq es lo mas seguro
      //sgg->thereAreMagneticMedia=.false.;
      //sgg->thereArePMLMagneticMedia=.false.;
      //medioespecial = .false.;
      //do ii=1,sgg->NumMedia
      //    buff='PEC media can only have index 0';
      //    if (sgg->Med(ii)->Is->PEC.or.sgg->Med(ii)->Is->Lumped) call STOPONERROR(layoutnumber,num_procs,buff);
      //    medioespecial =medioespecial .or. &
      //                   sgg->Med(ii)->Is->EDispersive   .or. &
      //                   sgg->Med(ii)->Is->multiport     .or. &
      //                   sgg->Med(ii)->Is->AnisMultiport .or. &
      //                   sgg->Med(ii)->Is->Anisotropic   .or. &
      //                   sgg->Med(ii)->Is->ThinSlot      .or. &
      //                   sgg->Med(ii)->Is->MDispersive;
      //end do
      //if (.not.medioespecial) then
      //    do ii=1,sgg->NumMedia
      //        if ((SGG->Med(ii)->mur >1.001).or.(SGG->Med(ii)->mur <0.999).or. &
      //            (abs(SGG->Med(ii)->sigmam) >1.0e-3_RKIND0)) then
      //            if (.not. sgg->Med(ii)->Is->PML ) then
      //                sgg->thereAreMagneticMedia=.true.;
      //            else
      //                sgg->thereArePMLMagneticMedia=.true.;
      //            end if
      //        end if
      //    end do
      //else
      sgg->thereAreMagneticMedia=.true.;
      sgg->thereArePMLMagneticMedia=.true.;
      //end if
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return;

      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine initConformalBoundingBox(sgg, bbox)
         type(SGGFDTDINFO_t), intent(in) :: sgg;
         type(XYZlimit_t), intent(inout) :: bbox;
         bbox->XI = -sgg->Alloc(iHx)->XI;
         bbox->XE = -sgg->Alloc(iHx)->XE;
         bbox->YI = -sgg->Alloc(iHy)->YI;
         bbox->YE = -sgg->Alloc(iHy)->YE;
         bbox->ZI = -sgg->Alloc(iHz)->ZI;
         bbox->ZE = -sgg->Alloc(iHz)->ZE;
      end subroutine

      function getDifferentEdgeRatios(conformal_media) result(res)
         type(ConformalMedia_t), dimension(:), allocatable, intent(in) :: conformal_media;
         integer :: i, j, k;
         real(kind=rkind), dimension(:), allocatable :: aux;
         real(kind=rkind), dimension(:), allocatable :: res;
         logical :: isNew;
         allocate(res(0));
         do

for (int i = 1; i <= conformal_media.extent(0); ++i) {
        for (int j = 1; j <= conformal_media(i).edge_media.extent(0); ++j) {
            bool isNew = true;
            for (int k = 1; k <= res.extent(0); ++k) {
                if (eq_ratio(res(k), conformal_media(i).edge_media(j).ratio, EDGE_RATIO_EQ_TOLERANCE)) {
                    isNew = false;
                }
            }
            if (isNew) {
                if (res.extent(0) == 0) {
                    res = array<real_kind>(1);
                    res(1) = conformal_media(i).edge_media(j).ratio;
                } else {
                    array<real_kind> aux(res.extent(0) + 1);
                    aux(1:res.extent(0)) = res;
                    aux(res.extent(0) + 1) = conformal_media(i).edge_media(j).ratio;
                    res = aux;
                }
            }
        }
    }
}

array<real_kind> getDifferentFaceRatios(const array<ConformalMedia_t, 1>& conformal_media) {
    array<real_kind> res(0);
    for (int i = 1; i <= conformal_media.extent(0); ++i) {
        for (int j = 1; j <= conformal_media(i).face_media.extent(0); ++j) {
            bool isNew = true;
            for (int k = 1; k <= res.extent(0); ++k) {
                if (eq_ratio(res(k), conformal_media(i).face_media(j).ratio, FACE_RATIO_EQ_TOLERANCE)) {
                    isNew = false;
                }
            }
            if (isNew) {
                if (res.extent(0) == 0) {
                    res = array<real_kind>(1);
                    res(1) = conformal_media(i).face_media(j).ratio;
                } else {
                    array<real_kind> aux(res.extent(0) + 1);
                    aux(1:res.extent(0)) = res;
                    aux(res.extent(0) + 1) = conformal_media(i).face_media(j).ratio;
                    res = aux;
                }
            }
        }
    }
    return res;
}

void addConformalMedia(SGGFDTDINFO_t& sgg, media_matrices_t& media, const ConformalMedia_t& conformal_media, 
                       const array<real_kind, 1>& edge_ratios, const array<real_kind, 1>& face_ratios, 
                       int contamedia, XYZlimit_t& bbox, const side_tris_map_t& side_map) {
    initConformalBoundingBox(sgg, bbox);

    array<real_kind, 1> edge_ratios_no_zero;
    if (findloc(edge_ratios, 0.0, 1) != 0) {
        edge_ratios_no_zero = array<real_kind>(edge_ratios.extent(0) - 1);
        int k = 0;
        for (int j = 1; j <= edge_ratios.extent(0); ++j) {
            if (edge_ratios(j) != 0) {
                k++;
                edge_ratios_no_zero(k) = edge_ratios(j);
            }
        }
    } else {
        edge_ratios_no_zero = edge_ratios;
    }

    array<real_kind, 1> face_ratios_no_zero;
    if (findloc(face_ratios, 0.0, 1) != 0) {
        face_ratios_no_zero = array<real_kind>(face_ratios.extent(0) - 1);
        int k = 0;
        for (int j = 1; j <= face_ratios.extent(0); ++j) {
            if (face_ratios(j) != 0) {
                k++;
                face_ratios_no_zero(k) = face_ratios(j);
            }
        }
    } else {
        face_ratios_no_zero = face_ratios;
    }

    int num_media = contamedia;
    addConformalEdgeMedia(sgg, media, conformal_media, num_media, edge_ratios_no_zero, bbox);
    num_media = contamedia + edge_ratios_no_zero.extent(0);
    addConformalFaceMedia(sgg, media, conformal_media, num_media, face_ratios_no_zero, bbox);
    addUndetectedBorderFaces(sgg, media, conformal_media, num_media, edge_ratios_no_zero, bbox, side_map);
}

void addConformalFaceMedia(SGGFDTDINFO_t& sgg, media_matrices_t& media, const ConformalMedia_t& conformal_media, 
                           int num_media, const array<real_kind, 1>& face_ratios, XYZlimit_t& bbox) {
    for (int j = 1; j <= conformal_media.n_faces_media; ++j) {
        int face_media = 0;
        if (conformal_media.face_media(j).ratio != 0) {
            face_media = num_media + findloc(face_ratios, conformal_media.face_media(j).ratio, 1);
            sgg.Med(face_media).Is.ConformalPEC = true;
            sgg.Med(face_media).Is.Needed = true;
            sgg.Med(face_media).Is.Volume = true;
            sgg.Med(face_media).Priority = prior_PEC;
            sgg.Med(face_media).Epr = this->mats.mats(1).eps / Eps0;
            sgg.Med(face_media).Sigma = 1.0e29_RKIND;
            sgg.Med(face_media).Mur = conformal_media.face_media(j).ratio * this->mats.mats(1).mu / Mu0;
            sgg.Med(face_media).SigmaM = 0.0_RKIND;
        }

        for (int k = 1; k <= conformal_media.face_media(j).n_elements; ++k) {
            array<int, 3> cell = conformal_media.face_media(j).faces(k).cell;
            if (cell(1) < bbox.xi) bbox.xi = cell(1);
            if (cell(1) > bbox.xe) bbox.xe = cell(1);
            if (cell(2) < bbox.yi) bbox.yi = cell(2);
            if (cell(2) > bbox.ye) bbox.ye = cell(2);
            if (cell(3) < bbox.zi) bbox.zi = cell(3);
            if (cell(3) > bbox.ze) bbox.ze = cell(3);

            switch (conformal_media.face_media(j).faces(k).direction) {
                case F_X:
                    media.sggMiHx(cell(1), cell(2), cell(3)) = face_media;
                    break;
                case F_Y:
                    media.sggMiHy(cell(1), cell(2), cell(3)) = face_media;
                    break;
                case F_Z:
                    media.sggMiHz(cell(1), cell(2), cell(3)) = face_media;
                    break;
            }
        }
    }
}

array<float, 3> getEdgeNormalFromTriangles(const array<triangle_t, 1>& triangles) {
    array<float, 3> res = {0.0f, 0.0f, 0.0f};
    for (int i = 1; i <= triangles.extent(0); ++i) {
        res = res + triangles(i).getNormal();
    }
    res = res / triangles.extent(0);
    return res;
}

void addConformalEdgeMedia(SGGFDTDINFO_t& sgg, media_matrices_t& media, const ConformalMedia_t& conformal_media, 
                           int num_media, const array<real_kind, 1>& edge_ratios, XYZlimit_t& bbox) {
    for (int j = 1; j <= conformal_media.n_edges_media; ++j) {
        int edge_media = 0;
        if (conformal_media.edge_media(j).ratio != 0) {
            edge_media = num_media + findloc(edge_ratios, conformal_media.edge_media(j).ratio, 1);
            sgg.Med(edge_media).Is.ConformalPEC = true;
            sgg.Med(edge_media).Is.Needed = true;
            sgg.Med(edge_media).Is.Volume = true;
            sgg.Med(edge_media).Priority = prior_PEC;
            sgg.Med(edge_media).Epr = (this->mats.mats(1).eps / conformal_media.edge_media(j).ratio) / Eps0;
            sgg.Med(edge_media).Sigma = 1.0e29_RKIND;
            sgg.Med(edge_media).Mur = this->mats.mats(1).mu / Mu0;
            sgg.Med(edge_media).SigmaM = 0.0_RKIND;
        }

        for (int k = 1; k <= conformal_media.edge_media(j).n_elements; ++k) {
            array<int, 3> cell = conformal_media.edge_media(j).edges(k).cell;

            if (cell(1) < bbox.xi) bbox.xi = cell(1);
            if (cell(1) > bbox.xe) bbox.xe = cell(1);
            if (cell(2) < bbox.yi) bbox.yi = cell(2);
            if (cell(2) > bbox.ye) bbox.ye = cell(2);
            if (cell(3) < bbox.zi) bbox.zi = cell(3);
            if (cell(3) > bbox.ze) bbox.ze = cell(3);

            switch (conformal_media.edge_media(j).edges(k).direction) {
                case E_X:
                    media.sggMiEx(cell(1), cell(2), cell(3)) = edge_media;
                    break;
                case E_Y:
                    media.sggMiEy(cell(1), cell(2), cell(3)) = edge_media;
                    break;
                case E_Z:
                    media.sggMiEz(cell(1), cell(2), cell(3)) = edge_media;
                    break;
            }
        }
    }
}

void addUndetectedBorderFaces(SGGFDTDINFO_t& sgg, media_matrices_t& media, const ConformalMedia_t& conformal_media, 
                              int num_media, const array<real_kind, 1>& edge_ratios, XYZlimit_t& bbox, 
                              const side_tris_map_t& side_map) {
    for (int j = 1; j <= conformal_media.n_edges_media; ++j) {
        if (conformal_media.edge_media(j).ratio == 0) {
            int edge_media = 0;
            for (int k = 1; k <= conformal_media.edge_media(j).n_elements; ++k) {
                array<int, 3> cell = conformal_media.edge_media(j).edges(k).cell;
                array<int, 4> key = {0, 0, 0, 0};
                key(1:3) = cell;

                switch (conformal_media.edge_media(j).edges(k).direction) {
                    case E_X:
                        key(4) = E_X;
                        if (side_map.hasKey(key)) {
                            array<triangle_t, 1> tris = side_map.getTrianglesFromSide(key);
                            array<float, 3> normal = getEdgeNormalFromTriangles(tris);
                            if (normal(2) < 0 && !sgg.med(media.sggMiHz(cell(1), cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHz(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2) + 1, cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1) + 1, cell(2), cell(3)) = edge_media;
                            } else if (normal(2) > 0 && !sgg.med(media.sggMiHz(cell(1), cell(2) - 1, cell(3))).is.conformalPEC) {
                                media.sggMiHz(cell(1), cell(2) - 1, cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2) - 1, cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2) - 1, cell(3)) = edge_media;
                                media.sggMiEy(cell(1) + 1, cell(2) - 1, cell(3)) = edge_media;
                            }

                            if (normal(3) < 0 && !sgg.med(media.sggMiHy(cell(1), cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHy(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3) + 1) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEz(cell(1) + 1, cell(2), cell(3)) = edge_media;
                            } else if (normal(3) > 0 && !sgg.med(media.sggMiHy(cell(1), cell(2), cell(3) - 1)).is.conformalPEC) {
                                media.sggMiHy(cell(1), cell(2), cell(3) - 1) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3) - 1) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3) - 1) = edge_media;
                                media.sggMiEz(cell(1) + 1, cell(2), cell(3) - 1) = edge_media;
                            }
                        }
                        break;
                    case E_Y:
                        key(4) = E_Y;
                        if (edge_media == 0 && side_map.hasKey(key)) {
                            array<triangle_t, 1> tris = side_map.getTrianglesFromSide(key);
                            array<float, 3> normal = getEdgeNormalFromTriangles(tris);
                            if (normal(3) < 0 && !sgg.med(media.sggMiHx(cell(1), cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3) + 1) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEz(cell(1), cell(2) + 1, cell(3)) = edge_media;
                            } else if (normal(3) > 0 && !sgg.med(media.sggMiHx(cell(1), cell(2), cell(3) - 1)).is.conformalPEC) {
                                media.sggMiHx(cell(1), cell(2), cell(3) - 1) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3) - 1) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3) - 1) = edge_media;
                                media.sggMiEz(cell(1), cell(2) + 1, cell(3) - 1) = edge_media;
                            }

                            if (normal(1) < 0 && !sgg.med(media.sggMiHz(cell(1), cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHz(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2) + 1, cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1) + 1, cell(2), cell(3)) = edge_media;
                            } else if (normal(1) > 0 && !sgg.med(media.sggMiHz(cell(1) - 1, cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHz(cell(1) - 1, cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1) - 1, cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1) - 1, cell(2) + 1, cell(3)) = edge_media;
                                media.sggMiEy(cell(1) - 1, cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3)) = edge_media;
                            }
                        }
                        break;
                    case E_Z:
                        key(4) = E_Z;
                        if (edge_media == 0 && side_map.hasKey(key)) {
                            array<triangle_t, 1> tris = side_map.getTrianglesFromSide(key);
                            array<float, 3> normal = getEdgeNormalFromTriangles(tris);
                            if (normal(2) < 0 && !sgg.med(media.sggMiHx(cell(1), cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2), cell(3) + 1) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEz(cell(1), cell(2) + 1, cell(3)) = edge_media;
                            } else if (normal(2) > 0 && !sgg.med(media.sggMiHx(cell(1), cell(2) - 1, cell(3))).is.conformalPEC) {
                                media.sggMiHx(cell(1), cell(2) - 1, cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2) - 1, cell(3)) = edge_media;
                                media.sggMiEy(cell(1), cell(2) - 1, cell(3) + 1) = edge_media;
                                media.sggMiEz(cell(1), cell(2) - 1, cell(3)) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3)) = edge_media;
                            }
                            if (normal(1) < 0 && !sgg.med(media.sggMiHy(cell(1), cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHy(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1), cell(2), cell(3) + 1) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3)) = edge_media;
                                media.sggMiEz(cell(1) + 1, cell(2), cell(3)) = edge_media;
                            } else if (normal(1) > 0 && !sgg.med(media.sggMiHy(cell(1) - 1, cell(2), cell(3))).is.conformalPEC) {
                                media.sggMiHy(cell(1) - 1, cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1) - 1, cell(2), cell(3)) = edge_media;
                                media.sggMiEx(cell(1) - 1, cell(2), cell(3) + 1) = edge_media;
                                media.sggMiEz(cell(1) - 1, cell(2), cell(3)) = edge_media;
                                media.sggMiEz(cell(1), cell(2), cell(3)) = edge_media;
                            }
                        }
                        break;
                }
            }
        }
    }
}

void read_TIMEFRECTRANSFsourcefiles(bool simu_devia) {
    real_kind deviafactor_Multiplier;
    real_kind unillo, tiempoant;
    if (simu_devia) {
        deviafactor_Multiplier = 1.0_RKIND;
    } else {
        deviafactor_Multiplier = 1.0_RKIND;
    }

    maxSourceValue = 0.0_RKIND;
    real_kind minSpaceStep = min(min(minval(sgg.dx), minval(sgg.dy)), minval(sgg.dz));

    for (int i = 1; i <= sgg.NumMedia; ++i) {
        if (sgg.Med(i).Is.ThinWire) {
            if (sgg.Med(i).wire(1).VsourceExists) {
                for (int CONTAVOLT = 1; CONTAVOLT <= sgg.Med(i).wire(1).NUMVOLTAGESOURCES; ++CONTAVOLT) {
                    bool errnofile = false;
                    bool exists = file_exists(trim(adjustl(sgg.Med(i).wire(1).VSource(CONTAVOLT).fichero.NAME)));
                    if (!exists) {
                        std::string buff = trim(adjustl(sgg.Med(i).wire(1).VSource(CONTAVOLT).fichero.name)) + " DOES NOT EXIST";
                        STOPONERROR(layoutnumber, num_procs, buff);
                    }
                    // Note: File opening and reading logic would require specific C++ file I/O implementation
                    // Placeholder for the logic:
                    // open(15, file=..., action='read')
                    // READ (15,*) tiempo1, field1
                    // READ (15,*) tiempo2, field2
                    // sgg.Med(i).wire(1).VSource(CONTAVOLT)...
                }
            }
        }
    }
}

VOLT).fichero.deltaSamples = tiempo2 - tiempo1;
            nsurfs = 3;
            //problemas con multivac
            // while (.not.eof(15))
            while (true) {
                READ (15,*, end=77) tiempo1, field1;
                if (field1/minspacestep > maxSourceValue) maxSourceValue=field1/minspacestep;
                nsurfs = nsurfs + 1;
            }
77:           continue;
            CLOSE (15);
            numus = nsurfs - 2;
            sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.NumSamples = numus;
            allocate(sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.Samples(0:numus));
            open(15, file=trim(adjustl(sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.NAME)),action='read');
            for (k = 0; k <= numus; k++) {
                tiempoant=tiempo1;
                READ (15,*) tiempo1, sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.Samples[k];
                sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.Samples[k] = sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.Samples[k] * 
                sgg->Med[i].wire[1].VSource[CONTAVOLT].Multiplier * 
                deviafactor_Multiplier;
                //evitar sampleos no uniformes
                if ((k>1)&&(k<numus-1)) {
                    unillo=(tiempo1-tiempoant)/sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.deltaSamples;
                    if ((unillo<0.9)||(unillo>1.0_RKIND/0.9)) {
                        if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) { //a tiempos muy altos ignoro el redondeo
                            buff=trim(adjustl(sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.NAME))+' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.';
                            if (!ignoresamplingerrors) {
                                CLOSE(15);
                                call STOPONERROR(layoutnumber,num_procs,buff);
                            }
                        }
                    }
                }
            }
            CLOSE (15);
            sgg->Med[i].wire[1].VSource[CONTAVOLT].fichero.NumSamples = numus;
        }
    }
}
if (sgg->Med[i].Is.ThinWire) {
    if (sgg->Med[i].wire[1].IsourceExists) {
        for (CONTACURR=1; CONTACURR<=sgg->Med[i].wire[1].NUMCURRENTSOURCES; CONTACURR++) {
            errnofile = false;
            inquire(file=trim(adjustl(sgg->Med[i].wire[1].ISource[CONTACURR].fichero.NAME)), EXIST=errnofile);
            if (!errnofile) {
                buff=trim(adjustl(sgg->Med[i].wire[1].ISource[CONTACURR].fichero.name))+' DOES NOT EXIST';
                call STOPONERROR(layoutnumber,num_procs,buff);
            }
            open(15, file=trim(adjustl(sgg->Med[i].wire[1].ISource[CONTACURR].fichero.NAME)),action='read');
            READ (15,*) tiempo1, field1;
            READ (15,*) tiempo2, field2;
            sgg->Med[i].wire[1].ISource[CONTACURR].fichero.deltaSamples = tiempo2 - tiempo1;
            nsurfs = 3;
            //problemas con multivac
            // while (.not.eof(15))
            while (true) {
                READ (15,*, end=79) tiempo1, field1;
                if (field1/minspacestep**2.0_RKIND > maxSourceValue) maxSourceValue=field1/minspacestep**2.0_RKIND; //aqui no tengo feeling, pero estas fuentes no se usan !!? repensar
                nsurfs = nsurfs + 1;
            }
79:           continue;
            CLOSE (15);
            numus = nsurfs - 2;
            sgg->Med[i].wire[1].ISource[CONTACURR].fichero.NumSamples = numus;
            allocate(sgg->Med[i].wire[1].ISource[CONTACURR].fichero.Samples(0:numus));
            open(15, file=trim(adjustl(sgg->Med[i].wire[1].ISource[CONTACURR].fichero.NAME)),action='read');
            for (k = 0; k <= numus; k++) {
                tiempoant=tiempo1;
                READ (15,*) tiempo1, sgg->Med[i].wire[1].ISource[CONTACURR].fichero.Samples[k];
                sgg->Med[i].wire[1].ISource[CONTACURR].fichero.Samples[k] = sgg->Med[i].wire[1].ISource[CONTACURR].fichero.Samples[k] * 
                sgg->Med[i].wire[1].ISource[CONTACURR].Multiplier * 
                deviafactor_Multiplier;
                //evitar sampleos no uniformes
                if ((k>1)&&(k<numus-1)) {
                    unillo=(tiempo1-tiempoant)/sgg->Med[i].wire[1].ISource[CONTACURR].fichero.deltaSamples;
                    if ((unillo<0.9)||(unillo>1.0_RKIND/0.9)) {
                        if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) { //a tiempos muy altos ignoro el redondeo
                            buff=trim(adjustl(sgg->Med[i].wire[1].ISource[CONTACURR].fichero.NAME))+' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.';
                            if (!ignoresamplingerrors) {
                                CLOSE(15);
                                call STOPONERROR(layoutnumber,num_procs,buff);
                            }
                        }
                    }
                }
            }
            CLOSE (15);
            sgg->Med[i].wire[1].ISource[CONTACURR].fichero.NumSamples = numus;
        }
    }
}
//
if (sgg->Med[i].Is.SlantedWire) {
    for (j = 1; j<=sgg->Med[i].SlantedWire[1].numNodes; j++) {
        if (sgg->Med[i].SlantedWire[1].nodes[j].VsourceExists) {
            errnofile = false;
            inquire(file=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.NAME)), EXIST=errnofile);
            if (!errnofile) {
                buff=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.name))+' DOES NOT EXIST';
                call STOPONERROR(layoutnumber,num_procs,buff);
            }
            open(15, file=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.NAME)),action='read');
            READ (15,*) tiempo1, field1;
            READ (15,*) tiempo2, field2;
            sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.deltaSamples = tiempo2 - tiempo1;
            nsurfs = 3;
            //problemas con multivac
            // while (.not.eof(15))
            while (true) {
                READ (15,*, end=179) tiempo1, field1;
                if (field1/minspacestep > maxSourceValue) maxSourceValue=field1/minspacestep;
                nsurfs = nsurfs + 1;
            }
179:          continue;
            CLOSE (15);
            numus = nsurfs - 2;
            sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.NumSamples = numus;
            allocate(sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.Samples(0:numus));
            open(15, file=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.NAME)),action='read');
            for (k = 0; k <= numus; k++) {
                tiempoant=tiempo1;
                READ (15,*) tiempo1, sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.Samples[k];
                sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.Samples[k] = sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.Samples[k] * 
                sgg->Med[i].SlantedWire[1].nodes[j].Vsource.Multiplier * 
                deviafactor_Multiplier;
                //evitar sampleos no uniformes
                if ((k>1)&&(k<numus-1)) {
                    unillo=(tiempo1-tiempoant)/sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.deltaSamples;
                    if ((unillo<0.9)||(unillo>1.0_RKIND/0.9)) {
                        if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) { //a tiempos muy altos ignoro el redondeo
                            buff=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.NAME))+' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.';
                            if (!ignoresamplingerrors) {
                                CLOSE(15);
                                call STOPONERROR(layoutnumber,num_procs,buff);
                            }
                        }
                    }
                }
            }
            CLOSE (15);
            sgg->Med[i].SlantedWire[1].nodes[j].Vsource.fichero.NumSamples = numus;
        }
        if (sgg->Med[i].SlantedWire[1].nodes[j].IsourceExists) {
            errnofile = false;
            inquire(file=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.NAME)), EXIST=errnofile);
            if (!errnofile) {
                buff=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.name))+' DOES NOT EXIST';
                call STOPONERROR(layoutnumber,num_procs,buff);
            }
            open(15, file=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.NAME)),action='read');
            READ (15,*) tiempo1, field1;
            READ (15,*) tiempo2, field2;
            sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.deltaSamples = tiempo2 - tiempo1;
            nsurfs = 3;
            //problemas con multivac
            // while (.not.eof(15))
            while (true) {
                READ (15,*, end=279) tiempo1, field1;
                if (field1/minspacestep**2.0_RKIND > maxSourceValue) maxSourceValue=field1/minspacestep**2.0_RKIND; //aqui no tengo feeling, pero estas fuentes no se usan !!? repensar
                nsurfs = nsurfs + 1;
            }
279:          continue;
            CLOSE (15);
            numus = nsurfs - 2;
            sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.NumSamples = numus;
            allocate(sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.Samples(0:numus));
            open(15, file=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.NAME)),action='read');
            for (k = 0; k <= numus; k++) {
                tiempoant=tiempo1;
                READ (15,*) tiempo1, sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.Samples[k];
                sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.Samples[k] = sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.Samples[k] * 
                sgg->Med[i].SlantedWire[1].nodes[j].Isource.Multiplier * 
                deviafactor_Multiplier;
                //evitar sampleos no uniformes
                if ((k>1)&&(k<numus-1)) {
                    unillo=(tiempo1-tiempoant)/sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.deltaSamples;
                    if ((unillo<0.9)||(unillo>1.0_RKIND/0.9)) {
                        if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) { //a tiempos muy altos ignoro el redondeo
                            buff=trim(adjustl(sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.NAME))+' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.';
                            if (!ignoresamplingerrors) {
                                CLOSE(15);
                                call STOPONERROR(layoutnumber,num_procs,buff);
                            }
                        }
                    }
                }
            }
            CLOSE (15);
            sgg->Med[i].SlantedWire[1].nodes[j].Isource.fichero.NumSamples = numus;
        }
    }
}
}
//nodal sources
for (j = 1; j <= sgg->NumNodalSources; j++) {
    //Read the time evolution file
    //
    if (!sgg->NodalSource[j].IsInitialValue) {
        errnofile = false;
        inquire(file=trim(adjustl(sgg->NodalSource[j].fichero.NAME)), EXIST=errnofile);
        if (!errnofile) {
            buff=trim(adjustl(sgg->NodalSource[j].fichero.name))+' DOES NOT EXIST';
            call STOPONERROR(layoutnumber,num_procs,buff);
        }
        open(15, file=trim(adjustl(sgg->NodalSource[j].fichero.NAME)),action='read');
        READ (15,*) tiempo1, field1;
        READ (15,*) tiempo2, field2;
        sgg->NodalSource[j].fichero.deltaSamples = tiempo2 - tiempo1;
        nsurfs = 3;
        while (true) {
            READ (15,*, end=78) tiempo1, field1;
            if (sgg->NodalSource[j].isElec) {
                if (field1 > maxSourceValue) maxSourceValue=field1;
            } else {
                if (field1*zvac > maxSourceValue) maxSourceValue=field1*zvac;
            }
            nsurfs = nsurfs + 1;
        }
78:           continue;
        CLOSE (15);
        numus = nsurfs - 2;
        allocate(sgg->NodalSource[j].fichero.Samples(0:numus));
        open(15, file=trim(adjustl(sgg->NodalSource[j].fichero.NAME)),action='read');
        for (k = 0; k <= numus; k++) {
            tiempoant=tiempo1;
            READ (15,*) tiempo1, sgg->NodalSource[j].fichero.Samples[k];
            sgg->NodalSource[j].fichero.Samples[k] = sgg->NodalSource[j].fichero.Samples[k] * deviafactor_Multiplier;
            //evitar sampleos no uniformes
            if ((k>1)&&(k<numus-1)) {
                unillo=(tiempo1-tiempoant)/sgg->NodalSource[j].fichero.deltaSamples;
                if ((unillo<0.9)||(unillo>1.0_RKIND/0.9)) {
                    if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) { //a tiempos muy altos ignoro el redondeo
                        buff=trim(adjustl(sgg->NodalSource[j].fichero.NAME))+' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.';
                        if (!ignoresamplingerrors) {
                            CLOSE(15);
                            call STOPONERROR(layoutnumber,num_procs,buff);
                        }
                    }
                }
            }
        }
        CLOSE (15);
        sgg->NodalSource[j].fichero.NumSamples = numus;
    } else { //es un initialvalue que no precisa fichero. incializar trivialmente
        numus = 0;
        sgg->NodalSource[j].fichero.deltaSamples = 1.0_RKIND; //must be non-null
        allocate(sgg->NodalSource[j].fichero.Samples(0:numus));
        for (k = 0; k <= numus; k++) {
            sgg->NodalSource[j].fichero.Samples[k] = 1.0_RKIND * deviafactor_Multiplier;
        }
        sgg->NodalSource[j].fichero.NumSamples = numus;
    }
}
//Plane wave sources
for (j = 1; j <= sgg->NumPlaneWaves; j++) {
    //Read the time evolution file
    //
    errnofile = false;
    inquire(file=trim(adjustl(sgg->PlaneWave[j].fichero.NAME)), EXIST=errnofile);
    if (!errnofile) {
        buff=trim(adjustl(sgg->PlaneWave[j].fichero.name))+' DOES NOT EXIST';
        call STOPONERROR(layoutnumber,num_procs,buff);
    }
    open(15, file=trim(adjustl(sgg->PlaneWave[j].fichero.NAME)),action='read');
    READ (15,*) tiempo1, field1;
    READ (15,*) tiempo2, field2;
    sgg->PlaneWave[j].fichero.deltaSamples = tiempo2 - tiempo1;
    nsurfs = 3;
    while (true) {
        READ (15,*, end=98) tiempo1, field1;
        if (field1 > maxSourceValue) maxSourceValue=field1;
        nsurfs = nsurfs + 1;
    }
98:           continue;
    CLOSE (15);
    numus = nsurfs - 2;
    allocate(sgg->PlaneWave[j].fichero.Samples(0:numus));
    open(15, file=trim(adjustl(sgg->PlaneWave[j].fichero.NAME)),action='read');
    for (k = 0; k <= numus; k++) {
        tiempoant=tiempo1;
        READ (15,*) tiempo1, sgg->PlaneWave[j].fichero.Samples[k];
        sgg->PlaneWave[j].fichero.Samples[k] = sgg->PlaneWave[j].fichero.Samples[k] * deviafactor_Multiplier;
        //evitar sampleos no uniformes
        if ((k>1)&&(k<numus-1)) {
            unillo=(tiempo1-tiempoant)/sgg->PlaneWave[j].fichero.deltaSamples;
            if ((unillo<0.9)||(unillo>1.0_RKIND/0.9)) {
                if (2.0_RKIND*(tiempo1-tiempoant)/(tiempo1+tiempoant)>1e-6_RKIND) { //a tiempos muy altos ignoro el redondeo
                    buff=trim(adjustl(sgg->PlaneWave[j].fichero.NAME))+' not uniformly sampled. Relaunch with -ignoresamplingerrors to ignore it.';
                    if (!ignoresamplingerrors) {
                        CLOSE(15);
                        call STOPONERROR(layoutnumber,num_procs,buff);
                    }
                }
            }
        }
    }
    CLOSE (15);
    sgg->PlaneWave[j].fichero.NumSamples = numus;
}
return;
}

!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!


void asignawiredisper(WireDispersiveParams_t& disp, const std::string& file) {

    int i, numPoles;
    double valreal, valimag;

    open (1712, file=file, action='READ');
    read (1712, *) numPoles;

    disp.NumPoles = numPoles;
    allocate(disp.res(1:numPoles));
    allocate(disp.p  (1:numPoles));
    for (i = 1; i <= numPoles; i++) {
        read (1712, *) valreal, valimag;
        disp.res[i] = CMPLX(valreal, valimag, CKIND);
    }
    for (i = 1; i <= numPoles; i++) {
        read (1712, *) valreal, valimag;
        disp.p[i] = CMPLX(valreal, valimag, CKIND);
    }
    read (1712, *) valreal, valimag;
    disp.d = CMPLX(valreal, valimag, CKIND);
    read (1712, *) valreal, valimag;
    disp.e = CMPLX(valreal, valimag, CKIND);

    return;
}

void asignadisper(FreqDepenMaterial_t* fdgeom) {

    if (fdgeom->l+fdgeom->LM !=0 ) {
        BUFF='ERROR: SECOND ORDER DISPERSIVE MEDIA UNSUPPORTED. TRANSLATE THEM TO FIRST ORDER ()';
        call WarnErrReport (buff,true);
    }
    //
    allocate(sgg->Med[contamedia].EDispersive(1));
    allocate(sgg->Med[contamedia].MDispersive(1));
    sgg->Med[contamedia].Priority = prior_FDB;
    sgg->Med[contamedia].Epr = fdgeom->eps11 / Eps0;
    sgg->Med[contamedia].Sigma = fdgeom->Sigma11;
    sgg->Med[contamedia].Mur = fdgeom->mu11 / Mu0;
    sgg->Med[contamedia].SigmaM = fdgeom->SigmaM11;
    //electric dispersion
    sgg->Med[contamedia].Edispersive[1].eps11=       fdgeom->e;

sgg.Med(contamedia).Edispersive(1).eps12 = fdgeom.eps12;
        sgg.Med(contamedia).Edispersive(1).eps13 = fdgeom.eps13;
        sgg.Med(contamedia).Edispersive(1).eps22 = fdgeom.eps22;
        sgg.Med(contamedia).Edispersive(1).eps23 = fdgeom.eps23;
        sgg.Med(contamedia).Edispersive(1).eps33 = fdgeom.eps33;

        sgg.Med(contamedia).Edispersive(1).mu11 = fdgeom.MU11;
        sgg.Med(contamedia).Edispersive(1).mu12 = fdgeom.MU12;
        sgg.Med(contamedia).Edispersive(1).mu13 = fdgeom.MU13;
        sgg.Med(contamedia).Edispersive(1).mu22 = fdgeom.MU22;
        sgg.Med(contamedia).Edispersive(1).mu23 = fdgeom.MU23;
        sgg.Med(contamedia).Edispersive(1).mu33 = fdgeom.MU33;

        sgg.Med(contamedia).EDispersive(1).SIGMA11 = fdgeom.SIGMA11;
        sgg.Med(contamedia).EDispersive(1).SIGMA12 = fdgeom.SIGMA12;
        sgg.Med(contamedia).EDispersive(1).SIGMA13 = fdgeom.SIGMA13;
        sgg.Med(contamedia).EDispersive(1).SIGMA22 = fdgeom.SIGMA22;
        sgg.Med(contamedia).EDispersive(1).SIGMA23 = fdgeom.SIGMA23;
        sgg.Med(contamedia).EDispersive(1).SIGMA33 = fdgeom.SIGMA33;

        sgg.Med(contamedia).EDispersive(1).SIGMAM11 = fdgeom.SIGMAM11;
        sgg.Med(contamedia).EDispersive(1).SIGMAM12 = fdgeom.SIGMAM12;
        sgg.Med(contamedia).EDispersive(1).SIGMAM13 = fdgeom.SIGMAM13;
        sgg.Med(contamedia).EDispersive(1).SIGMAM22 = fdgeom.SIGMAM22;
        sgg.Med(contamedia).EDispersive(1).SIGMAM23 = fdgeom.SIGMAM23;
        sgg.Med(contamedia).EDispersive(1).SIGMAM33 = fdgeom.SIGMAM33;

        sgg.Med(contamedia).Mdispersive(1).eps11 = sgg.Med(contamedia).Edispersive(1).eps11;
        sgg.Med(contamedia).Mdispersive(1).eps12 = sgg.Med(contamedia).Edispersive(1).eps12;
        sgg.Med(contamedia).Mdispersive(1).eps13 = sgg.Med(contamedia).Edispersive(1).eps13;
        sgg.Med(contamedia).Mdispersive(1).eps22 = sgg.Med(contamedia).Edispersive(1).eps22;
        sgg.Med(contamedia).Mdispersive(1).eps23 = sgg.Med(contamedia).Edispersive(1).eps23;
        sgg.Med(contamedia).Mdispersive(1).eps33 = sgg.Med(contamedia).Edispersive(1).eps33;

        sgg.Med(contamedia).Mdispersive(1).mu11 = sgg.Med(contamedia).Edispersive(1).mu11;
        sgg.Med(contamedia).Mdispersive(1).mu12 = sgg.Med(contamedia).Edispersive(1).mu12;
        sgg.Med(contamedia).Mdispersive(1).mu13 = sgg.Med(contamedia).Edispersive(1).mu13;
        sgg.Med(contamedia).Mdispersive(1).mu22 = sgg.Med(contamedia).Edispersive(1).mu22;
        sgg.Med(contamedia).Mdispersive(1).mu23 = sgg.Med(contamedia).Edispersive(1).mu23;
        sgg.Med(contamedia).Mdispersive(1).mu33 = sgg.Med(contamedia).Edispersive(1).mu33;

        sgg.Med(contamedia).MDISPERSIVE(1).SIGMA11 = sgg.Med(contamedia).EDispersive(1).SIGMA11;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMA12 = sgg.Med(contamedia).EDispersive(1).SIGMA12;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMA13 = sgg.Med(contamedia).EDispersive(1).SIGMA13;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMA22 = sgg.Med(contamedia).EDispersive(1).SIGMA22;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMA23 = sgg.Med(contamedia).EDispersive(1).SIGMA23;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMA33 = sgg.Med(contamedia).EDispersive(1).SIGMA33;

        sgg.Med(contamedia).MDISPERSIVE(1).SIGMAM11 = sgg.Med(contamedia).EDispersive(1).SIGMAM11;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMAM12 = sgg.Med(contamedia).EDispersive(1).SIGMAM12;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMAM13 = sgg.Med(contamedia).EDispersive(1).SIGMAM13;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMAM22 = sgg.Med(contamedia).EDispersive(1).SIGMAM22;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMAM23 = sgg.Med(contamedia).EDispersive(1).SIGMAM23;
        sgg.Med(contamedia).MDISPERSIVE(1).SIGMAM33 = sgg.Med(contamedia).EDispersive(1).SIGMAM33;

        // First order only. Second order terms do not play a role.
        sgg.Med(contamedia).EDispersive(1).NumPolRes11 = fdgeom.k11; // + fdgeom.l
        sgg.Med(contamedia).EDispersive(1).NumPolRes12 = fdgeom.k12;
        sgg.Med(contamedia).EDispersive(1).NumPolRes13 = fdgeom.k13;
        sgg.Med(contamedia).EDispersive(1).NumPolRes22 = fdgeom.k22;
        sgg.Med(contamedia).EDispersive(1).NumPolRes23 = fdgeom.k23;
        sgg.Med(contamedia).EDispersive(1).NumPolRes33 = fdgeom.k33;
        
        // Magnetic
        sgg.Med(contamedia).MDispersive(1).NumPolRes11 = fdgeom.KM11; // + fdgeom.LM
        sgg.Med(contamedia).MDispersive(1).NumPolRes12 = fdgeom.KM12;
        sgg.Med(contamedia).MDispersive(1).NumPolRes13 = fdgeom.KM13;
        sgg.Med(contamedia).MDispersive(1).NumPolRes22 = fdgeom.KM22;
        sgg.Med(contamedia).MDispersive(1).NumPolRes23 = fdgeom.KM23;
        sgg.Med(contamedia).MDispersive(1).NumPolRes33 = fdgeom.KM33;

        if (sgg.Med(contamedia).EDispersive(1).NumPolRes11 != 0) {
            sgg.Med(contamedia).Is.EDispersive = true;
            sgg.Med(contamedia).Is.EDispersiveANIS = false;
            sgg.Med(contamedia).Is.Dielectric = false;
        }
        if (sgg.Med(contamedia).EDispersive(1).NumPolRes12 + sgg.Med(contamedia).EDispersive(1).NumPolRes13 +
            sgg.Med(contamedia).EDispersive(1).NumPolRes22 + sgg.Med(contamedia).EDispersive(1).NumPolRes23 +
            sgg.Med(contamedia).EDispersive(1).NumPolRes33 != 0) {
            sgg.Med(contamedia).Is.EDispersive = true;
            sgg.Med(contamedia).Is.EDispersiveAnis = true;
            std::cout << "Error: anisotropic dispersive unsupported" << std::endl;
            std::exit(1);
            sgg.Med(contamedia).Is.Dielectric = false;
        }
        if (sgg.Med(contamedia).MDispersive(1).NumPolRes11 != 0) {
            sgg.Med(contamedia).Is.MDispersive = true;
            sgg.Med(contamedia).Is.Dielectric = false;
        }
        if (sgg.Med(contamedia).MDISPERSIVE(1).NumPolRes12 + sgg.Med(contamedia).MDISPERSIVE(1).NumPolRes13 +
            sgg.Med(contamedia).MDISPERSIVE(1).NumPolRes22 + sgg.Med(contamedia).MDISPERSIVE(1).NumPolRes23 +
            sgg.Med(contamedia).MDISPERSIVE(1).NumPolRes33 != 0) {
            sgg.Med(contamedia).Is.MDISPERSIVE = true;
            sgg.Med(contamedia).Is.MDISPERSIVEAnis = true;
            std::cout << "Error: anisotropic dispersive unsupported" << std::endl;
            std::exit(1);
            sgg.Med(contamedia).Is.Dielectric = false;
        }
        
        allocate(sgg.Med(contamedia).EDispersive(1).C11, sgg.Med(contamedia).EDispersive(1).a11,
                 sgg.Med(contamedia).EDispersive(1).NumPolRes11);
        allocate(sgg.Med(contamedia).EDispersive(1).C12, sgg.Med(contamedia).EDispersive(1).a12,
                 sgg.Med(contamedia).EDispersive(1).NumPolRes12);
        allocate(sgg.Med(contamedia).EDispersive(1).C13, sgg.Med(contamedia).EDispersive(1).a13,
                 sgg.Med(contamedia).EDispersive(1).NumPolRes13);
        allocate(sgg.Med(contamedia).EDispersive(1).C22, sgg.Med(contamedia).EDispersive(1).a22,
                 sgg.Med(contamedia).EDispersive(1).NumPolRes22);
        allocate(sgg.Med(contamedia).EDispersive(1).C23, sgg.Med(contamedia).EDispersive(1).a23,
                 sgg.Med(contamedia).EDispersive(1).NumPolRes23);
        allocate(sgg.Med(contamedia).EDispersive(1).C33, sgg.Med(contamedia).EDispersive(1).a33,
                 sgg.Med(contamedia).EDispersive(1).NumPolRes33);
        
        allocate(sgg.Med(contamedia).MDispersive(1).C11, sgg.Med(contamedia).MDispersive(1).a11,
                 sgg.Med(contamedia).MDispersive(1).NumPolRes11);
        allocate(sgg.Med(contamedia).MDispersive(1).C12, sgg.Med(contamedia).MDispersive(1).a12,
                 sgg.Med(contamedia).MDispersive(1).NumPolRes12);
        allocate(sgg.Med(contamedia).MDispersive(1).C13, sgg.Med(contamedia).MDispersive(1).a13,
                 sgg.Med(contamedia).MDispersive(1).NumPolRes13);
        allocate(sgg.Med(contamedia).MDispersive(1).C22, sgg.Med(contamedia).MDispersive(1).a22,
                 sgg.Med(contamedia).MDispersive(1).NumPolRes22);
        allocate(sgg.Med(contamedia).MDispersive(1).C23, sgg.Med(contamedia).MDispersive(1).a23,
                 sgg.Med(contamedia).MDispersive(1).NumPolRes23);
        allocate(sgg.Med(contamedia).MDispersive(1).C33, sgg.Med(contamedia).MDispersive(1).a33,
                 sgg.Med(contamedia).MDispersive(1).NumPolRes33);
        
        for (int k1 = 1; k1 <= fdgeom.k11; ++k1) {
            sgg.Med(contamedia).EDispersive(1).C11[k1] = fdgeom.a11[k1];
            sgg.Med(contamedia).EDispersive(1).a11[k1] = -fdgeom.b11[k1]; // The pole of ORIGINAL has a sign change
        }
        for (int k1 = 1; k1 <= fdgeom.k12; ++k1) {
            sgg.Med(contamedia).EDispersive(1).C12[k1] = fdgeom.a12[k1];
            sgg.Med(contamedia).EDispersive(1).a12[k1] = -fdgeom.b12[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.k13; ++k1) {
            sgg.Med(contamedia).EDispersive(1).C13[k1] = fdgeom.a13[k1];
            sgg.Med(contamedia).EDispersive(1).a13[k1] = -fdgeom.b13[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.k22; ++k1) {
            sgg.Med(contamedia).EDispersive(1).C22[k1] = fdgeom.a22[k1];
            sgg.Med(contamedia).EDispersive(1).a22[k1] = -fdgeom.b22[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.k23; ++k1) {
            sgg.Med(contamedia).EDispersive(1).C23[k1] = fdgeom.a23[k1];
            sgg.Med(contamedia).EDispersive(1).a23[k1] = -fdgeom.b23[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.k33; ++k1) {
            sgg.Med(contamedia).EDispersive(1).C33[k1] = fdgeom.a33[k1];
            sgg.Med(contamedia).EDispersive(1).a33[k1] = -fdgeom.b33[k1];
        }
        
        for (int k1 = 1; k1 <= fdgeom.KM11; ++k1) {
            sgg.Med(contamedia).MDispersive(1).C11[k1] = fdgeom.aM11[k1];
            sgg.Med(contamedia).MDispersive(1).a11[k1] = -fdgeom.bM11[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.KM12; ++k1) {
            sgg.Med(contamedia).MDispersive(1).C12[k1] = fdgeom.aM12[k1];
            sgg.Med(contamedia).MDispersive(1).a12[k1] = -fdgeom.bM12[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.KM13; ++k1) {
            sgg.Med(contamedia).MDispersive(1).C13[k1] = fdgeom.aM13[k1];
            sgg.Med(contamedia).MDispersive(1).a13[k1] = -fdgeom.bM13[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.KM22; ++k1) {
            sgg.Med(contamedia).MDispersive(1).C22[k1] = fdgeom.aM22[k1];
            sgg.Med(contamedia).MDispersive(1).a22[k1] = -fdgeom.bM22[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.KM23; ++k1) {
            sgg.Med(contamedia).MDispersive(1).C23[k1] = fdgeom.aM23[k1];
            sgg.Med(contamedia).MDispersive(1).a23[k1] = -fdgeom.bM23[k1];
        }
        for (int k1 = 1; k1 <= fdgeom.KM33; ++k1) {
            sgg.Med(contamedia).MDispersive(1).C33[k1] = fdgeom.aM33[k1];
            sgg.Med(contamedia).MDispersive(1).a33[k1] = -fdgeom.bM33[k1];
        }

        return;

    } // end subroutine asignadisper

    // ... (empty lines/comments preserved as structure markers if needed, but usually omitted in translation unless significant)

    } // end subroutine read_geomData

    // ... (empty lines/comments preserved as structure markers if needed)

    void read_limits_nogeom(int layoutnumber, int num_procs, SGGFDTDINFO_t& sgg, 
                            std::array<limit_t, 6>& fullsize, 
                            std::array<limit_t, 6>& SINPML_fullsize, 
                            const Parseador_t& this_obj, 
                            bool& MurAfterPML, bool& mur_exist) {
        // Note: In C++, we map the Fortran type structures. 
        // Assuming limit_t is a struct/class with appropriate members.
        // Assuming SGGFDTDINFO_t has the nested structure matching the Fortran derived types.
        // Assuming Parseador_t has 'despl' and 'general' and 'front' members.
        
        // dummy pointers from Fortran are not directly translatable without context of allocation strategy.
        // We assume sgg.dx, dy, dz are vectors or arrays allocated elsewhere or resized here.
        // Based on Fortran: allocate(sgg%dx(...)) implies dynamic allocation.
        
        // Displacement allocation
        // Fortran: allocate(sgg%dx(this%despl%mx1:this%despl%mx2-1), ...)
        // C++: Resize vectors to size (mx2 - mx1)
        int mx_size = this_obj.despl.mx2 - this_obj.despl.mx1;
        int my_size = this_obj.despl.my2 - this_obj.despl.my1;
        int mz_size = this_obj.despl.mz2 - this_obj.despl.mz1;

        sgg.dx.resize(mx_size);
        sgg.dy.resize(my_size);
        sgg.dz.resize(mz_size);

        // Material Matrix / Discretization X
        int tama = this_obj.despl.nx;
        if (tama == 1) {
            for (int i = 0; i < mx_size; ++i) { // Fortran index i from mx1 to mx2-1 maps to 0..mx_size-1 in vector
                sgg.dx[i] = this_obj.despl.desx[1]; // Fortran 1-based index
            }
        } else {
            if (tama != mx_size) {
                std::string buff = "Tamanio discretizacion distinto de la region";
                STOPONERROR(layoutnumber, num_procs, buff);
            }
            for (int i = 0; i < mx_size; ++i) {
                sgg.dx[i] = this_obj.despl.desx[i + 1]; // Fortran 1-based index
            }
        }

        // Discretization Y
        tama = this_obj.despl.nY;
        if (tama == 1) {
            for (int i = 0; i < my_size; ++i) {
                sgg.dy[i] = this_obj.despl.desY[1];
            }
        } else {
            if (tama != my_size) {
                std::string buff = "Tamanio discretizacion distinto de la region";
                STOPONERROR(layoutnumber, num_procs, buff);
            }
            for (int i = 0; i < my_size; ++i) {
                sgg.dy[i] = this_obj.despl.desY[i + 1];
            }
        }

        // Discretization Z
        tama = this_obj.despl.nZ;
        if (tama == 1) {
            for (int i = 0; i < mz_size; ++i) {
                sgg.dz[i] = this_obj.despl.desZ[1];
            }
        } else {
            if (tama != mz_size) {
                std::string buff = "Tamanio discretizacion distinto de la region";
                STOPONERROR(layoutnumber, num_procs, buff);
            }
            for (int i = 0; i < mz_size; ++i) {
                sgg.dz[i] = this_obj.despl.desZ[i + 1];
            }
        }

        // DISCRETIZATION LINES
        // X Lines
        tama = this_obj.despl.nx;
        std::vector<RKIND> lineasX(this_obj.despl.mx2 - this_obj.despl.mx1 + 1); // mx1 to mx2 inclusive
        // Fortran: lineasX(mx1) = originX + mx1 * desx(1)
        // C++: Index 0 corresponds to mx1
        lineasX[0] = this_obj.despl.originX + this_obj.despl.mx1 * this_obj.despl.desx[1];
        
        if (tama == 1) {
            for (int i = 0; i < mx_size; ++i) {
                // Fortran: lineasX(i+1) = desx(1) * (i-mx1+1) + lineasX(mx1)
                // C++: i is loop var 0..mx_size-1. Corresponds to Fortran i.
                // Index in vector: i+1.
                lineasX[i + 1] = this_obj.despl.desx[1] * (i - this_obj.despl.mx1 + 1) + lineasX[0];
            }
        } else {
            if (tama != mx_size) {
                std::string buff = "Tamanio discretizacion distinto de la region";
                STOPONERROR(layoutnumber, num_procs, buff);
            }
            for (int i = 0; i < mx_size; ++i) {
                // Fortran: lineasX(i+1) = desx(i) + lineasX(i)
                // C++: desx index i+1 (1-based)
                lineasX[i + 1] = this_obj.despl.desx[i + 1] + lineasX[i];
            }
        }

        // Y Lines
        tama = this_obj.despl.nY;
        std::vector<RKIND> lineasY(this_obj.despl.my2 - this_obj.despl.my1 + 1);
        lineasY[0] = this_obj.despl.originy + this_obj.despl.my1 * this_obj.despl.desY[1];
        
        if (tama == 1) {
            for (int i = 0; i < my_size; ++i) {
                lineasY[i + 1] = this_obj.despl.desY[1] * (i - this_obj.despl.my1 + 1) + lineasY[0];
            }
        } else {
            if (tama != my_size) {
                std::string buff = "Tamanio discretizacion distinto de la region";
                STOPONERROR(layoutnumber, num_procs, buff);
            }
            for (int i = 0; i < my_size; ++i) {
                lineasY[i + 1] = this_obj.despl.desY[i + 1] + lineasY[i];
            }
        }

        // Z Lines
        tama = this_obj.despl.nZ;
        std::vector<RKIND> lineasZ(this_obj.despl.mz2 - this_obj.despl.mz1 + 1);
        lineasZ[0] = this_obj.despl.originZ + this_obj.despl.mz1 * this_obj.despl.desZ[1];
        
        if (tama == 1) {
            for (int i = 0; i < mz_size; ++i) {
                lineasZ[i + 1] = this_obj.despl.desZ[1] * (i - this_obj.despl.mz1 + 1) + lineasZ[0];
            }
        } else {
            if (tama != mz_size) {
                std::string buff = "Tamanio discretizacion distinto de la region";
                STOPONERROR(layoutnumber, num_procs, buff);
            }
            for (int i = 0; i < mz_size; ++i) {
                lineasZ[i + 1] = this_obj.despl.desZ[i + 1] + lineasZ[i];
            }
        }

        // Allocate and assign LineX, LineY, LineZ
        // Fortran: allocate(sgg%LineX(mx1:mx2), ...)
        // C++: Resize vectors to size (mx2 - mx1 + 1)
        sgg.LineX.resize(this_obj.despl.mx2 - this_obj.despl.mx1 + 1);
        sgg.LineY.resize(this_obj.despl.my2 - this_obj.despl.my1 + 1);
        sgg.LineZ.resize(this_obj.despl.mz2 - this_obj.despl.mz1 + 1);

        // Fortran: sgg%LineX(mx1:mx2) = lineasX
        // C++: Copy entire vector
        sgg.LineX = lineasX;
        sgg.LineY = lineasY;
        sgg.LineZ = lineasZ;

        // General Parameters
        sgg.InitialTimeStep = 0;
        sgg.TimeSteps = this_obj.general.nmax;
        sgg.dt = this_obj.general.dt;

        // Border Initialization
        sgg.Border.IsBackPEC = false;
        sgg.Border.IsFrontPEC = false;
        sgg.Border.IsLeftPEC = false;
        sgg.Border.IsRightPEC = false;
        sgg.Border.IsUpPEC = false;
        sgg.Border.IsDownPEC = false;
        sgg.Border.IsBackPMC = false;
        sgg.Border.IsFrontPMC = false;
        sgg.Border.IsLeftPMC = false;
        sgg.Border.IsRightPMC = false;
        sgg.Border.IsUpPMC = false;
        sgg.Border.IsDownPMC = false;
        sgg.Border.IsBackPML = false;
        sgg.Border.IsFrontPML = false;
        sgg.Border.IsLeftPML = false;
        sgg.Border.IsRightPML = false;
        sgg.Border.IsUpPML = false;
        sgg.Border.IsDownPML = false;
        sgg.Border.IsBackPeriodic = false;
        sgg.Border.IsFrontPeriodic = false;
        sgg.Border.IsLeftPeriodic = false;
        sgg.Border.IsRightPeriodic = false;
        sgg.Border.IsUpPeriodic = false;
        sgg.Border.IsDownPeriodic = false;
        sgg.Border.IsBackMUR = false;
        sgg.Border.IsFrontMUR = false;
        sgg.Border.IsLeftMUR = false;
        sgg.Border.IsRightMUR = false;
        sgg.Border.IsUpMUR = false;
        sgg.Border.IsDownMUR = false;
        sgg.PML.NumLayers = 0;

        for (int i = 1; i <= 6; ++i) {
            if (this_obj.front.tipofrontera[i] == F_PML) {
                // Note: Fortran uses SELECT CASE on i. 
                // i=1: xmin (Back), i=2: xmax (Front), i=3: ymin (Left)
                // Assuming icoord, comi, jcoo, fine, etc. are defined in context or part of PML struct.
                // Since the code cuts off, I will translate up to the cut-off point.
                
                // The variables icoord, comi, jcoo, fine, etc. are likely local or member variables 
                // not shown in the snippet. I will assume they are accessible or part of the PML indexing logic.
                // Based on typical FDTD:
                // icoord/comi likely refer to X axis indices.
                // jcoo/fine likely refer to Y axis indices.
                
                // However, without their definition, I cannot complete the assignment correctly.
                // I will translate the structure up to the cut-off.
                
                if (i == 1) { // xmin
                    sgg.Border.IsBackPML = true;
                    // Assuming icoord and comi are defined in the surrounding scope or class
                    // sgg.PML.NumLayers(icoord, comi) = ...
                    // Since I don't have icoord/comi, I'll leave a placeholder or assume they are members.
                    // Given the previous context isn't here, I'll assume standard member access if they were members,
                    // but they look like local indices. 
                    // Let's assume they are part of the PML struct or global constants for this specific case.
                    // Actually, looking at Fortran, they are likely local integers defined earlier or members.
                    // I will assume they are members of sgg.PML or similar for the sake of syntax, 
                    // but strictly speaking, they are undefined here.
                    // I will translate the assignment assuming icoord/comi are available.
                    sgg.PML.NumLayers[/*icoord*/][/*comi*/] = this_obj.front.PROPIEDADESPML[i].NUMCAPAS;
                    sgg.PML.CoeffReflPML[/*icoord*/][/*comi*/] = this_obj.front.PROPIEDADESPML[i].REFL;
                    if (sgg.PML.CoeffReflPML[/*icoord*/][/*comi*/] >= 1.0_RKIND) {
                        sgg.PML.CoeffReflPML[/*icoord*/][/*comi*/] = 0.99999;
                    }
                    sgg.PML.orden[/*icoord*/][/*comi*/] = this_obj.front.PROPIEDADESPML[i].orden;
                } else if (i == 2) { // xmax
                    sgg.Border.IsFrontPML = true;
                    sgg.PML.NumLayers[/*icoord*/][/*fine*/] = this_obj.front.PROPIEDADESPML[i].NUMCAPAS;
                    sgg.PML.CoeffReflPML[/*icoord*/][/*fine*/] = this_obj.front.PROPIEDADESPML[i].REFL;
                    if (sgg.PML.CoeffReflPML[/*icoord*/][/*fine*/] >= 1.0_RKIND) {
                        sgg.PML.CoeffReflPML[/*icoord*/][/*fine*/] = 0.99999;
                    }
                    sgg.PML.orden[/*icoord*/][/*fine*/] = this_obj.front.PROPIEDADESPML[i].orden;
                } else if (i == 3) { // ymin
                    sgg.Border.IsLeftPML = true;
                    // Code cuts off here
                }
            }
        }
    }

rd, comi) = this->front->PROPIEDADESPML(i).NUMCAPAS;
               sgg->PML->CoeffReflPML(jcoord, comi) = this->front->PROPIEDADESPML(i).REFL;
               if (sgg->PML->CoeffReflPML(jcoord, comi) >= 1.0_RKIND) sgg->PML->CoeffReflPML(jcoord, comi) = 0.99999d0;
               sgg->PML->orden(jcoord, comi) = this->front->PROPIEDADESPML(i).orden;
               //ymax
             case (4):
               sgg->Border->IsRightPML = true;
               sgg->PML->NumLayers(jcoord, fine) = this->front->PROPIEDADESPML(i).NUMCAPAS;
               sgg->PML->CoeffReflPML(jcoord, fine) = this->front->PROPIEDADESPML(i).REFL;
               if (sgg->PML->CoeffReflPML(jcoord, fine) >= 1.0_RKIND) sgg->PML->CoeffReflPML(jcoord, fine) = 0.99999d0;
               sgg->PML->orden(jcoord, fine) = this->front->PROPIEDADESPML(i).orden;
               //zmin
             case (5):
               sgg->Border->IsDownPML = true;
               sgg->PML->NumLayers(kcoord, comi) = this->front->PROPIEDADESPML(i).NUMCAPAS;
               sgg->PML->CoeffReflPML(kcoord, comi) = this->front->PROPIEDADESPML(i).REFL;
               if (sgg->PML->CoeffReflPML(kcoord, comi) >= 1.0_RKIND) sgg->PML->CoeffReflPML(kcoord, comi) = 0.99999d0;
               sgg->PML->orden(kcoord, comi) = this->front->PROPIEDADESPML(i).orden;
               //zmax
             case (6):
               sgg->Border->IsUpPML = true;
               sgg->PML->NumLayers(kcoord, fine) = this->front->PROPIEDADESPML(i).NUMCAPAS;
               sgg->PML->CoeffReflPML(kcoord, fine) = this->front->PROPIEDADESPML(i).REFL;
               if (sgg->PML->CoeffReflPML(kcoord, fine) >= 1.0_RKIND) sgg->PML->CoeffReflPML(kcoord, fine) = 0.99999d0;
               sgg->PML->orden(kcoord, fine) = this->front->PROPIEDADESPML(i).orden;
            end switch
         } else if (this->front->tipofrontera(i) == F_MUR) {
            mur_exist = true;
            switch (i) {
               //xmin
             case (1):
               sgg->Border->IsBackMUR = true;
               //xmax
             case (2):
               sgg->Border->IsFrontMUR = true;
               //ymin
             case (3):
               sgg->Border->IsLeftMUR = true;
               //ymax
             case (4):
               sgg->Border->IsRightMUR = true;
               //zmin
             case (5):
               sgg->Border->IsDownMUR = true;
               //zmax
             case (6):
               sgg->Border->IsUpMUR = true;
            end switch
         } else if (this->front->tipofrontera(i) == F_PEC) {
            switch (i) {
               //xmin
             case (1):
               sgg->Border->IsBackPEC = true;
               //xmax
             case (2):
               sgg->Border->IsFrontPEC = true;
               //ymin
             case (3):
               sgg->Border->IsLeftPEC = true;
               //ymax
             case (4):
               sgg->Border->IsRightPEC = true;
               //zmin
             case (5):
               sgg->Border->IsDownPEC = true;
               //zmax
             case (6):
               sgg->Border->IsUpPEC = true;
            end switch
         } else if (this->front->tipofrontera(i) == F_PMC) {
            switch (i) {
               //xmin
             case (1):
               sgg->Border->IsBackPMC = true;
               //xmax
             case (2):
               sgg->Border->IsFrontPMC = true;
               //ymin
             case (3):
               sgg->Border->IsLeftPMC = true;
               //ymax
             case (4):
               sgg->Border->IsRightPMC = true;
               //zmin
             case (5):
               sgg->Border->IsDownPMC = true;
               //zmax
             case (6):
               sgg->Border->IsUpPMC = true;
            end switch
         } else if (this->front->tipofrontera(i) == F_Per) {
            switch (i) {
               //xmin
             case (1):
               sgg->Border->IsBackPeriodic = true;
               //xmax
             case (2):
               sgg->Border->IsFrontPeriodic = true;
               //ymin
             case (3):
               sgg->Border->IsLeftPeriodic = true;
               //ymax
             case (4):
               sgg->Border->IsRightPeriodic = true;
               //zmin
             case (5):
               sgg->Border->IsDownPeriodic = true;
               //zmax
             case (6):
               sgg->Border->IsUpPeriodic = true;
            end switch
         }
      }
      //assign limits
      for (field = iEx; field <= iHz; ++field) {
         SINPML_fullsize(field)->XI = this->despl->mx1;
         SINPML_fullsize(field)->YI = this->despl->my1;
         SINPML_fullsize(field)->ZI = this->despl->mz1;
         SINPML_fullsize(field)->XE = this->despl->mx2;
         SINPML_fullsize(field)->YE = this->despl->my2;
         SINPML_fullsize(field)->ZE = this->despl->mz2;
      }
      //adjust the endings
      SINPML_fullsize(iEx)->XE = SINPML_fullsize(iEx)->XE - 1;
      SINPML_fullsize(iEy)->YE = SINPML_fullsize(iEy)->YE - 1;
      SINPML_fullsize(iEz)->ZE = SINPML_fullsize(iEz)->ZE - 1;
      //
      //
      SINPML_fullsize(iHx)->YE = SINPML_fullsize(iHx)->YE - 1;
      SINPML_fullsize(iHx)->ZE = SINPML_fullsize(iHx)->ZE - 1;
      SINPML_fullsize(iHy)->ZE = SINPML_fullsize(iHy)->ZE - 1;
      SINPML_fullsize(iHy)->XE = SINPML_fullsize(iHy)->XE - 1;
      SINPML_fullsize(iHz)->XE = SINPML_fullsize(iHz)->XE - 1;
      SINPML_fullsize(iHz)->YE = SINPML_fullsize(iHz)->YE - 1;
      //
      for (field = iEx; field <= iHz; ++field) {
         fullsize(field)->XI = SINPML_fullsize(field)->XI - sgg->PML->NumLayers(icoord, comi);
         fullsize(field)->YI = SINPML_fullsize(field)->YI - sgg->PML->NumLayers(jcoord, comi);
         fullsize(field)->ZI = SINPML_fullsize(field)->ZI - sgg->PML->NumLayers(kcoord, comi);
         fullsize(field)->XE = SINPML_fullsize(field)->XE + sgg->PML->NumLayers(icoord, fine);
         fullsize(field)->YE = SINPML_fullsize(field)->YE + sgg->PML->NumLayers(jcoord, fine);
         fullsize(field)->ZE = SINPML_fullsize(field)->ZE + sgg->PML->NumLayers(kcoord, fine);
      }
      //
      //readjust mur boundaries if necessary

      sgg->Border->IsBackMUR = (sgg->Border->IsBackMUR) || (sgg->Border->IsBackPML && MurAfterPML);
      sgg->Border->IsFrontMUR = (sgg->Border->IsFrontMUR) || (sgg->Border->IsFrontPML && MurAfterPML);
      sgg->Border->IsLeftMUR = (sgg->Border->IsLeftMUR) || (sgg->Border->IsLeftPML && MurAfterPML);
      sgg->Border->IsRightMUR = (sgg->Border->IsRightMUR) || (sgg->Border->IsRightPML && MurAfterPML);
      sgg->Border->IsUpMUR = (sgg->Border->IsUpMUR) || (sgg->Border->IsUpPML && MurAfterPML);
      sgg->Border->IsDownMUR = (sgg->Border->IsDownMUR) || (sgg->Border->IsDownPML && MurAfterPML);


      //readjust space steps accordingly 140815 para que esten bien allocateados los dx
      // Discretization Lines Matrix Resizing to accomodate PML regions
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     allocate(DummyD(SINPML_fullsize(iHx)->XI:SINPML_fullsize(iHx)->XE-1));
      DummyD = sgg->dx;
      deallocate(sgg->dx);
      //

      sgg->allocDxI = -sgg->PML->NumLayers(1, 1) + SINPML_fullsize(iHx)->XI - 1 - 1;
      sgg->allocDxE = SINPML_fullsize(iHx)->XE + sgg->PML->NumLayers(1, 2) + 1 + 1;
     allocate(sgg->dx(sgg->allocDxI:sgg->allocDxE));
      //
      sgg->dx(SINPML_fullsize(iHx)->XI:SINPML_fullsize(iHx)->XE-1) = DummyD(SINPML_fullsize(iHx)->XI:SINPML_fullsize(iHx)->XE-1);
      deallocate(DummyD);
      //
     allocate(DummyD(SINPML_fullsize(iHy)->YI:SINPML_fullsize(iHy)->YE-1));
      DummyD = sgg->dy;
      deallocate(sgg->dy);
      //
      sgg->allocDyI = -sgg->PML->NumLayers(2, 1) + SINPML_fullsize(iHy)->YI - 1 - 1;
      sgg->allocDyE = SINPML_fullsize(iHy)->YE + sgg->PML->NumLayers(2, 2) + 1 + 1;
     allocate(sgg->dy(sgg->allocDyI:sgg->allocDyE));
      //
      sgg->dy(SINPML_fullsize(iHy)->YI:SINPML_fullsize(iHy)->YE-1) = DummyD(SINPML_fullsize(iHy)->YI:SINPML_fullsize(iHy)->YE-1);
      deallocate(DummyD);
      //
     allocate(DummyD(SINPML_fullsize(iHz)->ZI:SINPML_fullsize(iHz)->ZE-1));
      DummyD = sgg->dz;
      deallocate(sgg->dz);
      //
      sgg->allocDzI = -sgg->PML->NumLayers(3, 1) + SINPML_fullsize(iHz)->ZI - 1 - 1;
      sgg->allocDzE = SINPML_fullsize(iHz)->ZE + sgg->PML->NumLayers(3, 2) + 1 + 1;
     allocate(sgg->dz(sgg->allocDzI:sgg->allocDzE));
      sgg->dz(SINPML_fullsize(iHz)->ZI:SINPML_fullsize(iHz)->ZE-1) = DummyD(SINPML_fullsize(iHz)->ZI:SINPML_fullsize(iHz)->ZE-1);
      deallocate(DummyD);
      //
      delta = sgg->dx(SINPML_fullsize(iHx)->XI);
      for (i = SINPML_fullsize(iHx)->XI - 1; i >= SINPML_fullsize(iHx)->XI - 1 - sgg->PML->NumLayers(1, 1) - 1; i -= 1) {
         sgg->dx(i) = delta;
      }
      delta = sgg->dx(SINPML_fullsize(iHx)->XE-1);
      for (i = SINPML_fullsize(iHx)->XE; i <= SINPML_fullsize(iHx)->XE + sgg->PML->NumLayers(1, 2) + 1 + 1; ++i) {
         sgg->dx(i) = delta;
      }
      //
      delta = sgg->dy(SINPML_fullsize(iHy)->YI);
      for (j = SINPML_fullsize(iHy)->YI - 1; j >= SINPML_fullsize(iHy)->YI - 1 - sgg->PML->NumLayers(2, 1) - 1; j -= 1) {
         sgg->dy(j) = delta;
      }
      //
      delta = sgg->dy(SINPML_fullsize(iHy)->YE-1);
      for (j = SINPML_fullsize(iHy)->YE; j <= SINPML_fullsize(iHy)->YE + sgg->PML->NumLayers(2, 2) + 1 + 1; ++j) {
         sgg->dy(j) = delta;
      }
      //
      delta = sgg->dz(SINPML_fullsize(iHz)->ZI);
      for (k = SINPML_fullsize(iHz)->ZI - 1; k >= SINPML_fullsize(iHz)->ZI - 1 - sgg->PML->NumLayers(3, 1) - 1; k -= 1) {
         sgg->dz(k) = delta;
      }
      //
      delta = sgg->dz(SINPML_fullsize(iHz)->ZE-1);
      for (k = SINPML_fullsize(iHz)->ZE; k <= SINPML_fullsize(iHz)->ZE + sgg->PML->NumLayers(3, 2) + 1 + 1; ++k) {
         sgg->dz(k) = delta;
      }
      //DISCRETIZATION LINES (TO BE DEPRECATED IN A NEAR FUTURE, ONLY NEEDED BY THE PLANEWAVE CORNER ROUTINE)
      //
     allocate(DummyD(SINPML_fullsize(iHx)->XI:SINPML_fullsize(iHx)->XE));
      DummyD = sgg->LineX;
      deallocate(sgg->LineX);
     allocate(sgg->LineX(-sgg->PML->NumLayers(1, 1)+SINPML_fullsize(iHx)->XI-1:SINPML_fullsize(iHx)->XE+sgg->PML->NumLayers(1, 2)+1));
      sgg->LineX(SINPML_fullsize(iHx)->XI:SINPML_fullsize(iHx)->XE) = DummyD(SINPML_fullsize(iHx)->XI:SINPML_fullsize(iHx)->XE);
      deallocate(DummyD);
      //
     allocate(DummyD(SINPML_fullsize(iHy)->YI:SINPML_fullsize(iHy)->YE));
      DummyD = sgg->LineY;
      deallocate(sgg->LineY);
     allocate(sgg->LineY(-sgg->PML->NumLayers(2, 1)+SINPML_fullsize(iHy)->YI-1:SINPML_fullsize(iHy)->YE+sgg->PML->NumLayers(2, 2)+1));
      sgg->LineY(SINPML_fullsize(iHy)->YI:SINPML_fullsize(iHy)->YE) = DummyD(SINPML_fullsize(iHy)->YI:SINPML_fullsize(iHy)->YE);
      deallocate(DummyD);
      //
     allocate(DummyD(SINPML_fullsize(iHz)->ZI:SINPML_fullsize(iHz)->ZE));
      DummyD = sgg->LineZ;
      deallocate(sgg->LineZ);
     allocate(sgg->LineZ(-sgg->PML->NumLayers(3, 1)+SINPML_fullsize(iHz)->ZI-1:SINPML_fullsize(iHz)->ZE+sgg->PML->NumLayers(3, 2)+1));
      sgg->LineZ(SINPML_fullsize(iHz)->ZI:SINPML_fullsize(iHz)->ZE) = DummyD(SINPML_fullsize(iHz)->ZI:SINPML_fullsize(iHz)->ZE);
      deallocate(DummyD);
      //
      delta = sgg->LineX(SINPML_fullsize(iHx)->XI+1) - sgg->LineX(SINPML_fullsize(iHx)->XI);
      for (i = SINPML_fullsize(iHx)->XI - 1; i >= SINPML_fullsize(iHx)->XI - 1 - sgg->PML->NumLayers(1, 1); i -= 1) {
         sgg->LineX(i) = sgg->LineX(i+1) - delta;
      }
      delta = sgg->LineX(SINPML_fullsize(iHx)->XE) - sgg->LineX(SINPML_fullsize(iHx)->XE-1);
      for (i = SINPML_fullsize(iHx)->XE + 1; i <= SINPML_fullsize(iHx)->XE + sgg->PML->NumLayers(1, 2) + 1; ++i) {
         sgg->LineX(i) = sgg->LineX(i-1) + delta;
      }
      //
      delta = sgg->LineY(SINPML_fullsize(iHy)->YI+1) - sgg->LineY(SINPML_fullsize(iHy)->YI);
      for (j = SINPML_fullsize(iHy)->YI - 1; j >= SINPML_fullsize(iHy)->YI - 1 - sgg->PML->NumLayers(2, 1); j -= 1) {
         sgg->LineY(j) = sgg->LineY(j+1) - delta;
      }
      //
      delta = sgg->LineY(SINPML_fullsize(iHy)->YE) - sgg->LineY(SINPML_fullsize(iHy)->YE-1);
      for (j = SINPML_fullsize(iHy)->YE + 1; j <= SINPML_fullsize(iHy)->YE + sgg->PML->NumLayers(2, 2) + 1; ++j) {
         sgg->LineY(j) = sgg->LineY(j-1) + delta;
      }
      //
      delta = sgg->LineZ(SINPML_fullsize(iHz)->ZI+1) - sgg->LineZ(SINPML_fullsize(iHz)->ZI);
      for (k = SINPML_fullsize(iHz)->ZI - 1; k >= SINPML_fullsize(iHz)->ZI - 1 - sgg->PML->NumLayers(3, 1); k -= 1) {
         sgg->LineZ(k) = sgg->LineZ(k+1) - delta;
      }
      //
      delta = sgg->LineZ(SINPML_fullsize(iHz)->ZE) - sgg->LineZ(SINPML_fullsize(iHz)->ZE-1);
      for (k = SINPML_fullsize(iHz)->ZE + 1; k <= SINPML_fullsize(iHz)->ZE + sgg->PML->NumLayers(3, 2) + 1; ++k) {
         sgg->LineZ(k) = sgg->LineZ(k-1) + delta;
      }
      //2012
      //update actual number of media
      //!!!!!!

      //      do i = SINPML_fullsize(iHx)->XI - 1,  SINPML_fullsize(iHx)->XE + sgg->PML->NumLayers(1, 2) + 1
      //        write (6678,*) 'lineX ',i, sgg->LineX (i)
      //      end do
      //!!!
      //      do j = SINPML_fullsize(iHy)->YI - 1,  SINPML_fullsize(iHy)->YE + sgg->PML->NumLayers(2, 2) + 1
      //        write (6678,*) 'lineY ',j,sgg->LineY (j)
      //      end do
      //!!!
      //      do k = SINPML_fullsize(iHz)->ZI - 1, SINPML_fullsize(iHz)->ZE + sgg->PML->NumLayers(3, 2) + 1
      //        write (6678,*) 'lineZ ',k,sgg->LineZ (k)
      //      end do


      return;
      //
   }
   //



   //!!!!!!!!!!!!!!!PREPROCESADOR PARA SKIN-DEPTH 09/07/13
   void prepro_skindepth(this, fichin) {
      int pozi, tama, j, k;
      char multiportFile[BUFSIZE];
      Parseador_t& this_ref = *this;
      const char* fichin_ref = fichin;
      char restocadena[BUFSIZE];
      int my_iostat;

      open(7533, "UGRskindepthmatlab.layers");
      close(7533, status="delete");
      my_iostat = 0;
9306: if(my_iostat != 0) write(*, fmt="(a)", advance="no"), '.'; //if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9306,'.',quienmpi,'UGRskindepthmatlab.layers'
      open(7533, "UGRskindepthmatlab.layers", err=9306, iostat=my_iostat, status="new", action="write");

      tama = this->LossyThinSurfs->length;
      for (j = 1; j <= tama; ++j) {
         if (abs(this->LossyThinSurfs->cs(j)->SigmaM(1)) <= 1.0e-2_RKIND) { //SGBCs que hay que sustituir
            sprintf(multiportFile, "%s_z11.txt", trim(adjustl(this->LossyThinSurfs->cs(j)->files)));
            //
            //09/07/13 !los SGBCs con skindepth se deben preprocesar

            //crea el fichero de entrada para usar con el compilado de Matlab
            pozi = index(multiportFile, "_z11.txt");
            write(7533, "(a)") trim(adjustl(multiportFile[1:pozi-1]));
            write(7533, *) "layers", this->LossyThinSurfs->cs(j)->numcapas;
            for (k = 1; k <= this->LossyThinSurfs->cs(j)->numcapas; ++k) {
               write(7533, *) "eps", k, this->LossyThinSurfs->cs(j)->eps(k);
               write(7533, *) "mu", k, this->LossyThinSurfs->cs(j)->mu(k);
               write(7533, *) "sigma", k, this->LossyThinSurfs->cs(j)->sigma(k);
               write(7533, *) "thickness", k, this->LossyThinSurfs->cs(j)->thk(k);
            }
            write(7533, *) "fmin", 10**4;
            write(7533, *) "fmax", 10**9;
            write(7533, *) "order", 24;
         }
      }
      close(7533);
      return;

   } //prepro_skindepth

   void AssigLossyOrPECtoNodes(sgg, media) {
      SGGFDTDINFO_t& sgg_ref = *sgg;
      media_matrices_t& media_ref = *media;

      bool ispec, isSGBC, IsComposite, islossy, input_conformal_flag, NODALMENTEIGUALES, iguaSGM, iguaSIG, iguaMUR, iguaPEC, iguaLOS, iguaEPR, ISconformal;
      double sigt, epst, SIGMA, SIGMAM, EPR, MUR;
      int i, j, k, n, kmenos1, jmenos1, imenos1, med[6], r, imed, i1;
      char buff[BUFSIZE];

      for (k = sgg_ref.Alloc(iEz)->ZI; k <= sgg_ref.Alloc(iEz)->ZE; ++k) {
         for (j = sgg_ref.Alloc(iEy)->YI; j <= sgg_ref.Alloc(iEy)->YE; ++j) {
            for (i = sgg_ref.Alloc(iEx)->XI; i <= sgg_ref.Alloc(iEx)->XE; ++i) {
               imenos1 = i - 1;
               jmenos1 = j - 1;
               kmenos1 = k - 1;
               if (i - 1 < sgg_ref.alloc(iEx)->XI) imenos1 = i;
               if (j - 1 < sgg_ref.alloc(iEy)->YI) jmenos1 = j;
               if (k - 1 < sgg_ref.alloc(iEz)->ZI) kmenos1 = k;

               med[0] = media_ref.sggMiEx(i, j, k);
               med[1] = media_ref.sggMiEx(imenos1, j, k);
               med[2] = media_ref.sggMiEy(i, j, k);
               med[3] = media_ref.sggMiEy(i, jmenos1, k);
               med[4] = media_ref.sggMiEz(i, j, k);
               med[5] = media_ref.sggMiEz(i, j, kmenos1);
               sigma = 0.0_RKIND;
               sigmam = 0.0_RKIND;
               epr = 0.0_RKIND;
               mur = 0.0_RKIND;
               ISPEC = false;
               ISLOSSY = false;
               for (i1 = 0; i1 <= 5; ++i1) {
                  imed = med[i1];
                  sigma = max(sigma, sgg_ref.Med(imed)->sigma);
                  sigmam = max(sigmam, sgg_ref.Med(imed)->sigmam);
                  epr = epr + sgg_ref.Med(imed)->epr / 6.0_RKIND;
                  mur = mur + sgg_ref.Med(imed)->mur / 6.0_RKIND;
                  if ((sgg_ref.med(imed)->is->PEC) || (imed == 0)) isPEC = true;
               }
               if ((!isPEC) && (sigma >= 1e-4)) {
                  islossy = true;
               } else {
                  islossy = false;
               }
               // CREAR NUEVO MEDIO y asignarle sus propiedades de acuerdo a sus adyacencias
               if (!(sgg_ref.med(med[0])->is->PML || sgg_ref.med(med[1])->is->PML || sgg_ref.med(med[2])->is->PML || sgg_ref.med(med[3])->is->PML || sgg_ref.med(med[4])->is->PML || sgg_ref.med(med[5])->is->PML)) {
                  if ((MED(0) != MED(1)) || (MED(1) != MED(2)) || (MED(2) != MED(3)) || (MED(3) != MED(4)) || (MED(4) != MED(5)) || (MED(5) != MED(0))) {
                     NODALMENTEIGUALES = false;
                     busqueda: for (I1 = 0; I1 <= SGG->NUMMEDIA; ++I1) {
                        //cambios 230817 bug milano borja en rutina iguales
                        iguaSGM = IGUALES(SGG->MED(I1)->SIGMAM, SIGMAM);
                        iguaSIG = IGUALES(SGG->MED(I1)->SIGMA, SIGMA); //sgg230817 al poner sigma 1e29 y pec ademas, la rutina de iguales fallaba
                        iguaEPR = IGUALES(SGG->MED(I1)->EPR, EPR);
                        iguaMUR = IGUALES(SGG->MED(I1)->MUR, MUR);
                        iguaPEC = (SGG->MED(I1)->iS->PEC == ISPEC);
                        iguaLOS = (SGG->MED(I1)->iS->LOSSY == ISLOSSY);
                        ISconformal = ((SGG->MED(I1)->iS->already_YEEadvanced_byconformal) || (SGG->MED(I1)->is->split_and_useless));
                        NODALMENTEIGUALES = NODALMENTEIGUALES || (((iguaSGM && iguaSIG && iguaEPR && iguaMUR && iguaLOS) || iguaPEC) && (!ISconformal));
                        if (nodalmenteiguales) break busqueda;
                     } end loop busqueda
                     if (!NODALMENTEIGUALES) {
                        if (SGG->NUMMEDIA + 1 > SGG->ALLOCmed) {
                           READJUST(SGG->ALLOCmed, sgg->med, 2 * SGG->ALLOCmed); //LO HAgo REallocatando al doble. gENERO NUEVO PARAMETRO sgg%ALLOCmed. Pero esto es un guirigay.... 261115
                        }
                        SGG->NUMMEDIA = SGG->NUMMEDIA + 1;
                        media_ref.sggMiNo(i, j, k) = SGG->NUMMEDIA;
                        r = media_ref.sggMiNo(i, j, k);
                        sgg->med(r)->sigma = sigma;
                        sgg->med(r)->sigmam = sigmam;
                        sgg->med(r)->epr = epr;
                        sgg->med(r)->mur = mur;
                        sgg->med(r)->is->PEC = ISPEC; //ojo con estos medios que el sistema de prioridades ya no les afecta porque esta rutina va despues del preprocess 03116
                        sgg->med(r)->is->LOSSY = ISLOSSY;
                        sgg->med(r)->is->needed = true; //sgg 220817 por defecto lo he puesto en readjust a false
                        //write(113,*) '.NOT.NODALMENTEIGUALES--> ',i,j,k,' - ',med(0),med(1),med(2),med(3),med(4),med(5),' - ',SGG->NUMMEDIA
                     } else {
                        media_ref.sggMiNo(i, j, k) = i1; //PUEDE QUE NO SEAN IGUALES PERO NODALMENTE LO SON (SOLO A E

bool IGUALES(real A, real B) {
    real ERR;
    bool igual = false;
    if (std::abs(A + B) > 1e-20) {
        ERR = 2.0_RKIND * std::abs((A - B) / (A + B));
        if (ERR < 1e-2_RKIND) igual = true; // en tanto por ciento me apanio con un 1 por ciento
    } else {
        ERR = std::abs(A - B);
        if (ERR < 1e-20_RKIND) igual = true; // en valor absoluto para valores casi nulos le pido que el error sea casi nulo
    }
    return igual;
}

void populatePlaneWaveRC(planeonde_t& PlaneWave, bool simu_devia) {
    int kkk;
    real theta, phi, alpha, beta, alpha1, alpha2, amplitud, FACTOR;
    complex beta1, beta2;
    bool primeravez;
    char buff[BUFSIZE];

    call_random_seed();

    amplitud = 1.0;
    for (kkk = 1; kkk <= PlaneWave.numModes; kkk++) {
        primeravez = true;
        theta = 0.; phi = 0.;
        while (((2.0_RKIND * pi * std::sin(theta) < phi)) || primeravez) { // moglie
            primeravez = false;
            call_random_number(theta);
            theta = pi * theta;
            call_random_number(phi);
            phi = 2.0_RKIND * pi * phi;
        } // moglie
        phi = phi / std::sin(theta); // moglie

        // ahora la polarizacion
        // si los hago asi hay apegotonamiento en los polos 281115 pero con los beta tampoco me sale. Seguir pensando y hacerlo con poincare algun dia 281115
        // generado con ortogonalidad_teM_parafuentesRC.nb
        // ojo que el atan de fortran y de mathematica estan invertidos!!!!
        
        call_random_number(beta);
        beta = 2.0_RKIND * pi * beta;
        
        alpha1 = atan2(
            Cos(theta) / Sqrt(Cos(theta)**2.0_RKIND + Cos(beta - phi)**2.0 * Sin(theta)**2.0),
            -((Cos(beta - phi) * Sin(theta)) / Sqrt(Cos(theta)**2.0_RKIND + Cos(beta - phi)**2.0 * Sin(theta)**2.0))
        );
        alpha2 = atan2(
            -(Cos(theta) / Sqrt(Cos(theta)**2.0_RKIND + Cos(beta - phi)**2.0 * Sin(theta)**2.0)),
            (Cos(beta - phi) * Sin(theta)) / Sqrt(Cos(theta)**2.0_RKIND + Cos(beta - phi)**2.0 * Sin(theta)**2.0)
        );
        
        if ((alpha1 <= pi) && (alpha1 >= 0.0)) {
            alpha = alpha1;
        } else if ((alpha2 <= pi) && (alpha2 >= 0.0)) {
            alpha = alpha2;
        } else {
            goto label_2;
            // write(buff,*) 'Error generando direcciones aleatorias para RC planewaves. '
            // call STOPONERROR(0,0,buff)
        }
        
label_2:
        // ahora la incertumbre en la posicion
        call_random_number(FACTOR);
        PlaneWave.INCERT[kkk] = PlaneWave.incertMax * FACTOR;

        PlaneWave.px[kkk] = Sin(theta) * Cos(phi);
        PlaneWave.py[kkk] = Sin(theta) * Sin(phi);
        PlaneWave.pz[kkk] = Cos(theta);
        PlaneWave.ex[kkk] = amplitud * Sin(alpha) * Cos(beta);
        PlaneWave.ey[kkk] = amplitud * Sin(alpha) * Sin(beta);
        PlaneWave.ez[kkk] = amplitud * Cos(alpha);

        if (std::abs(PlaneWave.px[kkk]**2. + PlaneWave.py[kkk]**2. + PlaneWave.pz[kkk]**2. - 1.) > 1e-4) {
            goto label_1;
        }
        if (std::abs(PlaneWave.ex[kkk]**2. + PlaneWave.ey[kkk]**2. + PlaneWave.ez[kkk]**2. - amplitud**2.) > 1e-4) {
            goto label_1;
        }
        if (std::abs(PlaneWave.px[kkk] * PlaneWave.ex[kkk] + PlaneWave.py[kkk] * PlaneWave.ey[kkk] + PlaneWave.pz[kkk] * PlaneWave.ez[kkk]) >= 1e-4) {
            goto label_1;
        }
    }

label_1:;

    int thefileno = 888;
    open(thefileno, "rc_EP.dat", FORM='formatted');
    while (true) {
        int status = read(thefileno, "(i5,12e19.9e3)", kkk, PlaneWave.px[kkk], PlaneWave.py[kkk], PlaneWave.pz[kkk],
            PlaneWave.ex[kkk], PlaneWave.ey[kkk], PlaneWave.ez[kkk], PlaneWave.INCERT[kkk]);
        if (status != 0) break;
    }
    close(thefileno);

    if (!simu_devia) { // solo lo escribe el principal
        thefileno = 888;
        open(thefileno, "rc_EP.dat", FORM='formatted');
        for (kkk = 1; kkk <= PlaneWave.numModes; kkk++) {
            write(thefileno, "(i5,12e19.9e3)") kkk, PlaneWave.px[kkk], PlaneWave.py[kkk], PlaneWave.pz[kkk],
                PlaneWave.ex[kkk], PlaneWave.ey[kkk], PlaneWave.ez[kkk], PlaneWave.INCERT[kkk];
        }
        close(thefileno);
    }
}

void cuentatags(Parseador_t& this, tagtype_t tagtype, int layoutnumber, const char* fichin) {
    bool foundDuplicate;
    int numertag, i, j, k, m, tama, tama2, tama3, tama2p, tama3p, precounting, acum, thefileno;
    char tagToCheck[BUFSIZE];

    for (precounting = 0; precounting <= 1; precounting++) {
        numertag = 0;
        
        tama = this.pecregs.nvols;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            if (i > 1) {
                if (this.pecregs.vols[i].tag == this.pecregs.vols[i-1].tag) { // do not increase
                    numertag = numertag - 1;
                }
            }
            if (precounting == 1) tagtype.tag[numertag] = this.pecregs.vols[i].tag;
        }
        
        tama = this.pecregs.nsurfs;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            if (i > 1) {
                if (this.pecregs.surfs[i].tag == this.pecregs.surfs[i-1].tag) { // do not increase
                    numertag = numertag - 1;
                }
            }
            if (precounting == 1) tagtype.tag[numertag] = this.pecregs.surfs[i].tag;
        }
        
        tama = this.pecregs.nLINS;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            if (i > 1) {
                if (this.pecregs.lins[i].tag == this.pecregs.lins[i-1].tag) { // do not increase
                    numertag = numertag - 1;
                }
            }
            if (precounting == 1) tagtype.tag[numertag] = this.pecregs.lins[i].tag;
        }

        tama = this.pmcregs.nvols;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            if (i > 1) {
                if (this.pmcregs.vols[i].tag == this.pmcregs.vols[i-1].tag) { // do not increase
                    numertag = numertag - 1;
                }
            }
            if (precounting == 1) tagtype.tag[numertag] = this.pmcregs.vols[i].tag;
        }
        
        tama = this.pmcregs.nsurfs;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            if (i > 1) {
                if (this.pmcregs.surfs[i].tag == this.pmcregs.surfs[i-1].tag) { // do not increase
                    numertag = numertag - 1;
                }
            }
            if (precounting == 1) tagtype.tag[numertag] = this.pmcregs.surfs[i].tag;
        }
        
        tama = this.pmcregs.nLINS;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            if (i > 1) {
                if (this.pmcregs.lins[i].tag == this.pmcregs.lins[i-1].tag) { // do not increase
                    numertag = numertag - 1;
                }
            }
            if (precounting == 1) tagtype.tag[numertag] = this.pmcregs.lins[i].tag;
        }

        tama = this.DielRegs.nvols;
        for (i = 1; i <= tama; i++) {
            tama2 = this.DielRegs.vols[i].n_c1P;
            tama3 = this.DielRegs.vols[i].n_c2P;

            if ((tama2 == 0) && (tama3 == 0)) {
                print *,"Bug in Dielectric Tags. Missing coordinates";
                stop;
            }

            // Check c1P coordinates
            checkDielectricComponentTags(this.DielRegs.vols[i], this.DielRegs.vols(1:i-1), i-1, 
                "c1P", numertag, tagtype, precounting, 
                "Bug in Dielectric Volume Tags");

            // Check c2P coordinates
            checkDielectricComponentTags(this.DielRegs.vols[i], this.DielRegs.vols(1:i-1), i-1, 
                "c2P", numertag, tagtype, precounting, 
                "Bug in Dielectric Volume Tags");
        }

        // Similar for surfaces
        tama = this.DielRegs.nSurfs;
        for (i = 1; i <= tama; i++) {
            tama2 = this.DielRegs.surfs[i].n_c1P;
            tama3 = this.DielRegs.surfs[i].n_c2P;

            if ((tama2 == 0) && (tama3 == 0)) {
                print *,"Bug in Dielectric Tags. Missing coordinates";
                stop;
            }

            // Check c1P coordinates
            checkDielectricComponentTags(this.DielRegs.surfs[i], this.DielRegs.surfs(1:i-1), i-1, 
                "c1P", numertag, tagtype, precounting, 
                "Bug in Dielectric Surface Tags");

            // Check c2P coordinates
            checkDielectricComponentTags(this.DielRegs.surfs[i], this.DielRegs.surfs(1:i-1), i-1, 
                "c2P", numertag, tagtype, precounting, 
                "Bug in Dielectric Surface Tags");
        }

        // Similar for surfaces (lines)
        tama = this.DielRegs.nLins;
        for (i = 1; i <= tama; i++) {
            tama2 = this.DielRegs.Lins[i].n_c1P;
            tama3 = this.DielRegs.Lins[i].n_c2P;

            if ((tama2 == 0) && (tama3 == 0)) {
                print *,"Bug in Dielectric Tags. Missing coordinates";
                stop;
            }

            // Check c1P coordinates
            checkDielectricComponentTags(this.DielRegs.lins[i], this.DielRegs.lins(1:i-1), i-1, 
                "c2P", numertag, tagtype, precounting, 
                "Bug in Dielectric Surface Tags");

            // Check c2P coordinates
            checkDielectricComponentTags(this.DielRegs.lins[i], this.DielRegs.lins(1:i-1), i-1, 
                "c2P", numertag, tagtype, precounting, 
                "Bug in Dielectric Surface Tags");
        }

        tama = this.animats.nvols;
        for (i = 1; i <= tama; i++) {
            tama2 = this.DielRegs.vols[i].n_c1P;
            tama3 = this.DielRegs.vols[i].n_c2P;

            if ((tama2 == 0) && (tama3 == 0)) {
                print *,"Bug in Animat Tags. Missing coordinates";
                stop;
            }

            // Check c1P coordinates
            checkAnimatedComponentTags(this.aniMats.vols[i], this.aniMats.vols(1:i-1), i-1, 
                "c1P", numertag, tagtype, precounting, 
                "Bug in Animat Volume Tags");

            // Check c2P coordinates
            checkAnimatedComponentTags(this.aniMats.vols[i], this.aniMats.vols(1:i-1), i-1, 
                "c2P", numertag, tagtype, precounting, 
                "Bug in Animat Volume Tags");
        }
        
        tama = this.animats.nSurfs;
        for (i = 1; i <= tama; i++) {
            tama2 = this.DielRegs.Surfs[i].n_c1P;
            tama3 = this.DielRegs.Surfs[i].n_c2P;

            if ((tama2 == 0) && (tama3 == 0)) {
                print *,"Bug in Animat Tags. Missing coordinates";
                stop;
            }

            // Check c1P coordinates
            checkAnimatedComponentTags(this.aniMats.Surfs[i], this.aniMats.Surfs(1:i-1), i-1, 
                "c1P", numertag, tagtype, precounting, 
                "Bug in Animat Surface Tags");

            // Check c2P coordinates
            checkAnimatedComponentTags(this.aniMats.Surfs[i], this.aniMats.Surfs(1:i-1), i-1, 
                "c2P", numertag, tagtype, precounting, 
                "Bug in Animat Surface Tags");
        }
        
        tama = this.animats.nLins;
        for (i = 1; i <= tama; i++) {
            tama2 = this.DielRegs.Lins[i].n_c1P;
            tama3 = this.DielRegs.Lins[i].n_c2P;

            if ((tama2 == 0) && (tama3 == 0)) {
                print *,"Bug in Animat Tags. Missing coordinates";
                stop;
            }

            // Check c1P coordinates
            checkAnimatedComponentTags(this.aniMats.Lins[i], this.aniMats.Lins(1:i-1), i-1, 
                "c1P", numertag, tagtype, precounting, 
                "Bug in Animat Line Tags");

            // Check c2P coordinates
            checkAnimatedComponentTags(this.aniMats.Lins[i], this.aniMats.Lins(1:i-1), i-1, 
                "c2P", numertag, tagtype, precounting, 
                "Bug in Animat Line Tags");
        }

        tama = this.frqdepmats.nvols;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            tama2 = this.frqdepmats.vols[i].n_c;
            if (tama2 != 0) {
                if (i > 1) {
                    if (this.frqdepmats.vols[i].c(1).tag == this.frqdepmats.vols[i-1].c(1).tag) { // do not increase
                        numertag = numertag - 1;
                    }
                }
            }
            if (precounting == 1) {
                if (tama2 != 0) {
                    tagtype.tag[numertag] = this.frqdepmats.vols[i].c(1).tag;
                } else {
                    print *,"bug in tags. ";
                    stop;
                }
                for (j = 1; j <= tama2; j++) {
                    if (trim(adjustl(this.frqdepmats.vols[i].c(j).tag)) != trim(adjustl(tagtype.tag[numertag]))) {
                        print *,"bug in tags. ";
                        stop;
                    }
                }
            }
        }
        
        tama = this.frqdepmats.nsurfs;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            tama2 = this.frqdepmats.surfs[i].n_c;
            if (tama2 != 0) {
                if (i > 1) {
                    if (this.frqdepmats.surfs[i].c(1).tag == this.frqdepmats.surfs[i-1].c(1).tag) { // do not increase
                        numertag = numertag - 1;
                    }
                }
            }
            if (precounting == 1) {
                if (tama2 != 0) {
                    tagtype.tag[numertag] = this.frqdepmats.surfs[i].c(1).tag;
                } else {
                    print *,"bug in tags. ";
                    stop;
                }
                for (j = 1; j <= tama2; j++) {
                    if (trim(adjustl(this.frqdepmats.surfs[i].c(j).tag)) != trim(adjustl(tagtype.tag[numertag]))) {
                        print *,"bug in tags. ";
                        stop;
                    }
                }
            }
        }
        
        tama = this.frqdepmats.nlins;
        for (i = 1; i <= tama; i++) {
            numertag = numertag + 1;
            tama2 = this.frqdepmats.lins[i].n_c;
            if (tama2 != 0) {
                if (i > 1) {
                    if (this.frqdepmats.lins[i].c(1).tag == this.frqdepmats.lins[i-1].c(1).tag) { // do not increase
                        numertag = numertag - 1;
                    }
                }
            }
            if (precounting == 1) {
                if (tama2 != 0) {
                    tagtype.tag[numertag] = this.LossyThinSurfs.cs(i).C(1).tag;
                } else {
                    print *,"bug in tags. ";
                    stop;
                }
                for (j = 1; j <= tama2; j++) {
                    if (trim(adjustl(this.frqdepmats.lins[i].c(j).tag)) != trim(adjustl(tagtype.tag[numertag]))) {
                        print *,"bug in tags. ";
                        stop;
                    }
                }
            }
        }
    }
}

}
      }
      //!!!
      tama = this->LossyThinSurfs.length;
      for (i = 1; i <= tama; ++i) {
         checkLossyTags(this->LossyThinSurfs.cs[i],
                        this->LossyThinSurfs.cs,
                        i - 1, numertag, tagtype, precounting);
      }

      tama = this->twires.n_tw;
      for (i = 1; i <= tama; ++i) {
         numertag = numertag + 1;
         tama2 = this->twires.TW[i].N_TWC;
         if (tama2 != 0) {
            if (i > 1) {
               if (this->twires.TW[i].TWC[1].tag == this->twires.TW[i - 1].TWC[1].tag) { // do not increase
                  numertag = numertag - 1;
               }
            }
         }
         if (precounting == 1) {
            if (tama2 != 0) {
               tagtype.tag[numertag] = this->twires.TW[i].TWC[1].tag;
            } else {
               std::cout << "bug in tags. " << std::endl;
               std::exit(1);
            }
            for (j = 1; j <= tama2; ++j) {
               if (trim(adjustl(this->twires.TW[i].TWC[j].tag)) != trim(adjustl(tagtype.tag[numertag]))) {
                  std::cout << "bug in tags. " << std::endl;
                  std::exit(1);
               }
            }
         }
      }
#ifdef CompileWithMTLN
      {
         cable_t* ptr;
         for (i = 1; i <= this->mtln.n_unsh + this->mtln.n_sh; ++i) {
            ptr = this->mtln.cables[i].ptr;
            if (auto* unshielded_ptr = dynamic_cast<unshielded_multiwire_t*>(ptr)) {
               numertag = numertag + 1;
               if (precounting == 1) {
                  tagtype.tag[numertag] = this->mtln.cables[i].ptr->tag;
               }
            }
         }
      }
#endif

      tama = this->swires.n_sw;
      for (i = 1; i <= tama; ++i) {
         numertag = numertag + 1;
         tama2 = this->swires.SW[i].n_swc;
         if (tama2 != 0) {
            if (i > 1) {
               if (this->swires.SW[i].swc[1].tag == this->swires.SW[i - 1].swc[1].tag) { // do not increase
                  numertag = numertag - 1;
               }
            }
         }
         if (precounting == 1) {
            if (tama2 != 0) {
               tagtype.tag[numertag] = this->swires.SW[i].swc[1].tag;
            } else {
               std::cout << "bug in tags. " << std::endl;
               std::exit(1);
            }
            for (j = 1; j <= tama2; ++j) {
               if (trim(adjustl(this->swires.SW[i].swc[j].tag)) != trim(adjustl(tagtype.tag[numertag]))) {
                  std::cout << "bug in tags. " << std::endl;
                  std::exit(1);
               }
            }
         }
      }

      tama = this->tSlots.n_tg;
      for (i = 1; i <= tama; ++i) {
         numertag = numertag + 1;
         tama2 = this->tSlots.Tg[i].N_tgc;
         if (tama2 != 0) {
            if (i > 1) {
               if (this->tSlots.Tg[i].TgC[1].tag == this->tSlots.Tg[i - 1].TgC[1].tag) { // do not increase
                  numertag = numertag - 1;
               }
            }
         }
         if (precounting == 1) {
            if (tama2 != 0) {
               tagtype.tag[numertag] = this->tSlots.Tg[i].TgC[1].tag;
            } else {
               std::cout << "bug in tags. " << std::endl;
               std::exit(1);
            }
            for (j = 1; j <= tama2; ++j) {
               if (trim(adjustl(this->tSlots.Tg[i].TgC[j].tag)) != trim(adjustl(tagtype.tag[numertag]))) {
                  std::cout << "bug in tags. " << std::endl;
                  std::exit(1);
               }
            }
         }
      }

      if (this->conformalRegs.volumes != nullptr) {
         numertag = numertag + 1;
         if (precounting == 1) tagtype.tag[numertag] = this->conformalRegs.volumes[1].tag;
      }
      if (this->conformalRegs.surfaces != nullptr) {
         numertag = numertag + 1;
         if (precounting == 1) tagtype.tag[numertag] = this->conformalRegs.surfaces[1].tag;
      }
      //!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!
      if (precounting == 0) {
         tagtype.numertags = numertag;
         tagtype.tag.resize(numertag + 1); // uno mas para luego jugar
         for (size_t k = 0; k < tagtype.tag.size(); ++k) tagtype.tag[k] = "";
      } else { // elimina repetidos
         for (i = 1; i <= numertag; ++i) {
            for (j = i + 1; j <= numertag; ++j) {
               if (trim(adjustl(tagtype.tag[i])) == trim(adjustl(tagtype.tag[j]))) {
                  tagtype.tag[j] = "";
               }
            }
         }
         i = 1;
         acum = 0;
         while ((i <= numertag) && (acum <= numertag + 1)) {
            if (trim(adjustl(tagtype.tag[i])) == "") {
               // Shift array left
               for (int k = i; k < numertag; ++k) {
                  tagtype.tag[k] = tagtype.tag[k + 1];
               }
               tagtype.tag[numertag] = ""; // Clear last element to avoid garbage
               acum = acum + 1;
            } else {
               i = i + 1;
            }
         }
         numertag = i - 1;
         tagtype.numertags = numertag;
      }
   } // do !del precounting


   return;
}

void checkDielectricComponentTags(Dielectric_t& component, const std::vector<Dielectric_t>& prev_components, int n_prev, const std::string& coord_type, int& numertag, tagtype_t& tagtype, int precounting, const std::string& error_msg) {
   int tama2;

   if (coord_type == "c1P") {
      tama2 = component.n_c1P;
   } else {
      tama2 = component.n_c2P;
   }

   for (int j = 1; j <= tama2; ++j) {
      numertag = numertag + 1;
      checkDielectricTagForDuplicate(component, prev_components, n_prev, j, coord_type, numertag, tagtype, precounting, error_msg);
   }
}

void checkDielectricTagForDuplicate(Dielectric_t& component, const std::vector<Dielectric_t>& prev_components, int n_prev, int idx, const std::string& coord_type, int& numertag, tagtype_t& tagtype, int precounting, const std::string& error_msg) {
   bool foundDuplicate = false;
   std::string tagToCheck;
   int k, m;

   if (coord_type == "c1P") {
      tagToCheck = trim(adjustl(component.c1P[idx].tag));
   } else {
      tagToCheck = trim(adjustl(component.c2P[idx].tag));
   }

   if (idx > 1) {
      for (k = 1; k <= idx - 1; ++k) {
         if (coord_type == "c1P") {
            if (tagToCheck == trim(adjustl(component.c1P[k].tag))) {
               foundDuplicate = true;
               break;
            }
         } else {
            if (tagToCheck == trim(adjustl(component.c2P[k].tag))) {
               foundDuplicate = true;
               break;
            }
         }
      }
   }

   if (!foundDuplicate) {
      for (m = 1; m <= n_prev; ++m) {
         if (prev_components[m].n_c1P > 0) {
            for (k = 1; k <= prev_components[m].n_c1P; ++k) {
               if (tagToCheck == trim(adjustl(prev_components[m].c1P[k].tag))) {
                  std::cout << error_msg << std::endl;
                  std::cout << "Duplicate tag found: " << tagToCheck << std::endl;
                  std::exit(1);
               }
            }
         }

         if (prev_components[m].n_c2P > 0) {
            for (k = 1; k <= prev_components[m].n_c2P; ++k) {
               if (tagToCheck == trim(adjustl(prev_components[m].c2P[k].tag))) {
                  std::cout << error_msg << std::endl;
                  std::cout << "Duplicate tag found: " << tagToCheck << std::endl;
                  std::exit(1);
               }
            }
         }
      }
   }

   if (foundDuplicate) {
      numertag = numertag - 1;
   } else if (precounting == 1) {
      tagtype.tag[numertag] = tagToCheck;
   }
}

void checkAnimatedComponentTags(ANISOTROPICbody_t& component, const std::vector<ANISOTROPICbody_t>& prev_components, int n_prev, const std::string& coord_type, int& numertag, tagtype_t& tagtype, int precounting, const std::string& error_msg) {
   int tama2;

   if (coord_type == "c1P") {
      tama2 = component.n_c1P;
   } else {
      tama2 = component.n_c2P;
   }

   for (int j = 1; j <= tama2; ++j) {
      numertag = numertag + 1;
      checkAnimatedTagForDuplicate(component, prev_components, n_prev, j, coord_type, numertag, tagtype, precounting, error_msg);
   }
}

void checkAnimatedTagForDuplicate(ANISOTROPICbody_t& component, const std::vector<ANISOTROPICbody_t>& prev_components, int n_prev, int idx, const std::string& coord_type, int& numertag, tagtype_t& tagtype, int precounting, const std::string& error_msg) {
   bool foundDuplicate = false;
   std::string tagToCheck;
   int k, m;

   if (coord_type == "c1P") {
      tagToCheck = trim(adjustl(component.c1P[idx].tag));
   } else {
      tagToCheck = trim(adjustl(component.c2P[idx].tag));
   }

   if (idx > 1) {
      for (k = 1; k <= idx - 1; ++k) {
         if (coord_type == "c1P") {
            if (tagToCheck == trim(adjustl(component.c1P[k].tag))) {
               foundDuplicate = true;
               break;
            }
         } else {
            if (tagToCheck == trim(adjustl(component.c2P[k].tag))) {
               foundDuplicate = true;
               break;
            }
         }
      }
   }

   if (!foundDuplicate) {
      for (m = 1; m <= n_prev; ++m) {
         if (prev_components[m].n_c1P > 0) {
            for (k = 1; k <= prev_components[m].n_c1P; ++k) {
               if (tagToCheck == trim(adjustl(prev_components[m].c1P[k].tag))) {
                  std::cout << error_msg << std::endl;
                  std::cout << "Duplicate tag found: " << tagToCheck << std::endl;
                  std::exit(1);
               }
            }
         }

         if (prev_components[m].n_c2P > 0) {
            for (k = 1; k <= prev_components[m].n_c2P; ++k) {
               if (tagToCheck == trim(adjustl(prev_components[m].c2P[k].tag))) {
                  std::cout << error_msg << std::endl;
                  std::cout << "Duplicate tag found: " << tagToCheck << std::endl;
                  std::exit(1);
               }
            }
         }
      }
   }

   if (foundDuplicate) {
      numertag = numertag - 1;
   } else if (precounting == 1) {
      tagtype.tag[numertag] = tagToCheck;
   }
}

void checkLossyTags(LossyThinSurface_t& component, const std::vector<LossyThinSurface_t>& prev_components, int n_prev, int& numertag, tagtype_t& tagtype, int precounting) {
   int j;

   if (component.nc == 0) {
      std::cout << "Bug in LossyThinSurf Tags. Missing coordinates" << std::endl;
      std::exit(1);
   }

   for (j = 1; j <= component.nc; ++j) {
      numertag = numertag + 1;
      checkLossyTagForDuplicate(component, prev_components, n_prev, j, numertag, tagtype, precounting);
   }
}

void checkLossyTagForDuplicate(LossyThinSurface_t& component, const std::vector<LossyThinSurface_t>& prev_components, int n_prev, int idx, int& numertag, tagtype_t& tagtype, int precounting) {
   bool foundDuplicate = false;
   std::string tagToCheck;
   int k, m;

   tagToCheck = trim(adjustl(component.C[idx].tag));

   if (idx > 1) {
      for (k = 1; k <= idx - 1; ++k) {
         if (tagToCheck == trim(adjustl(component.C[k].tag))) {
            foundDuplicate = true;
            break;
         }
      }
   }

   if ((!foundDuplicate) && (n_prev > 0)) {
      for (m = 1; m <= n_prev; ++m) {
         if (prev_components[m].nc > 0) {
            for (k = 1; k <= prev_components[m].nc; ++k) {
               if (tagToCheck == trim(adjustl(prev_components[m].C[k].tag))) {
                  foundDuplicate = true;
               }
            }
         }
      }
   }

   if (foundDuplicate) {
      numertag = numertag - 1;
   } else if (precounting == 1) {
      tagtype.tag[numertag] = tagToCheck;
   }
}

int searchtag(const tagtype_t& tagtype, const std::string& tag) {
   int numertag = -1;
   for (int i = 1; i <= tagtype.numertags; ++i) {
      if (trim(adjustl(tagtype.tag[i])) == trim(adjustl(tag))) {
         numertag = i;
         break;
      }
   }
   return numertag;
}