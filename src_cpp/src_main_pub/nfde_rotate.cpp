```cpp
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>

// Forward declarations for types defined in NFDETypes_m
// In a real scenario, these would be in a header file "NFDETypes_m.h"
struct Parseador_t;
struct Desplazamiento_t;
struct MatrizMedios_t;
struct mtln_t;
struct FronteraPML_t;
struct PlaneWaves_t;
struct Boxes_t;
struct BoxesRegs_t;
struct DielRegs_t;
struct ThinWires_t;
struct SlantedWiresInfo_t;
struct ThinSlots_t;
struct LossyThinSurfs_t;
struct FRQDEPMATS_t;
struct FarField_Sonda_t;
struct Electric_Sonda_t;
struct Magnetic_Sonda_t;
struct coords_t;
struct BloqueProbe_t;
struct FreqDepenMaterial_t;
struct coords_scaled_t;

// Assuming these constants/enums are defined in NFDETypes_m
// They are used in the code but their definitions are not provided in the snippet.
// We declare them here as extern or assume they are available.
extern int iEx;
extern int iEy;
extern int iEz;
extern int NP_COR_EX;
extern int NP_COR_EY;
extern int NP_COR_EZ;
extern int NP_COR_hX;
extern int NP_COR_hY;
extern int NP_COR_hZ;
extern int iEXc;
extern int iEYc;
extern int iEZc;
extern int ihXc;
extern int ihYc;
extern int ihZc;
extern int iCurX;
extern int iCurY;
extern int iCurZ;
extern int iMEC;
extern int iMHC;
extern int iCur;

// RK and RKIND are likely kind parameters for real numbers
// Assuming RK is double and RKIND is a suffix or constant for double precision
using RK = double;
constexpr RK RKIND = 1.0; // Placeholder, usually handled by compiler kind selectors

namespace nfde_rotate_m {

    // Helper function for sign, mimicking Fortran sign(x, y)
    inline int sign(int x, int y) {
        return (y >= 0) ? std::abs(x) : -std::abs(x);
    }

    inline void nfde_rotate(Parseador_t& this_obj, int mpidir) {
        rotate_generateSpaceSteps(this_obj, mpidir);
        rotate_generateCurrent_Field_Sources(this_obj, mpidir);
        rotate_generatePlaneWaves(this_obj, mpidir);
        rotate_generateBoxSources(this_obj, mpidir);
        rotate_generateFronteras(this_obj, mpidir);
        rotate_generatePECs(this_obj, mpidir);
        rotate_generatePMCs(this_obj, mpidir);
        rotate_generateNONMetals(this_obj, mpidir);
        rotate_generateANISOTROPICs(this_obj, mpidir);
        rotate_generateThinWires(this_obj, mpidir);
        rotate_generateSlantedWires(this_obj, mpidir);
        rotate_generateThinSlots(this_obj, mpidir);
        rotate_generateLossyThinSurface(this_obj, mpidir);
        rotate_generateFDMs(this_obj, mpidir);
        rotate_generateSONDAs(this_obj, mpidir);
        rotate_generateMasSondas(this_obj, mpidir);
        rotate_generateBloqueProbes(this_obj, mpidir);
        rotate_generateVolumicProbes(this_obj, mpidir);
#ifdef CompileWithMTLN
        rotate_mtln(this_obj, mpidir);
#endif
    }

#ifdef CompileWithMTLN
    void rotate_mtln(Parseador_t& this_obj, int mpidir) {
        mtln_t* old_mtln = new mtln_t(this_obj.mtln); // Assuming copy constructor exists or implementing deep copy logic if needed
        // Note: Fortran 'allocate(..., source=...)' performs a copy.
        
        int i, j;
        int x, y, z, or;

        for (i = 0; i < static_cast<int>(old_mtln->cables.size()); ++i) {
            for (j = 0; j < static_cast<int>(old_mtln->cables[i].ptr->segments.size()); ++j) {
                x = old_mtln->cables[i].ptr->segments[j].x;
                y = old_mtln->cables[i].ptr->segments[j].y;
                z = old_mtln->cables[i].ptr->segments[j].z;
                or = old_mtln->cables[i].ptr->segments[j].orientation;

                if (mpidir == 2) {
                    this_obj.mtln.cables[i].ptr->segments[j].x = z;
                    this_obj.mtln.cables[i].ptr->segments[j].y = x;
                    this_obj.mtln.cables[i].ptr->segments[j].z = y;
                    switch (std::abs(or)) {
                        case 1:
                            this_obj.mtln.cables[i].ptr->segments[j].orientation = sign(2, or);
                            break;
                        case 2:
                            this_obj.mtln.cables[i].ptr->segments[j].orientation = sign(3, or);
                            break;
                        case 3:
                            this_obj.mtln.cables[i].ptr->segments[j].orientation = sign(1, or);
                            break;
                    }
                } else if (mpidir == 1) {
                    this_obj.mtln.cables[i].ptr->segments[j].x = y;
                    this_obj.mtln.cables[i].ptr->segments[j].y = z;
                    this_obj.mtln.cables[i].ptr->segments[j].z = x;
                    switch (std::abs(old_mtln->cables[i].ptr->segments[j].orientation)) {
                        case 1:
                            this_obj.mtln.cables[i].ptr->segments[j].orientation = sign(3, or);
                            break;
                        case 2:
                            this_obj.mtln.cables[i].ptr->segments[j].orientation = sign(1, or);
                            break;
                        case 3:
                            this_obj.mtln.cables[i].ptr->segments[j].orientation = sign(2, or);
                            break;
                    }
                }
            }
        }
        delete old_mtln;
    }
#endif

    void rotate_generateSpaceSteps(Parseador_t& this_obj, int mpidir) {
        Desplazamiento_t* old_despl = new Desplazamiento_t(this_obj.despl);
        MatrizMedios_t* old_matriz = new MatrizMedios_t(this_obj.matriz);

        int oxi, oyi, ozi;
        double roxi, royi, rozi;
        std::vector<double>* poxi = nullptr;
        std::vector<double>* poyi = nullptr;
        std::vector<double>* pozi = nullptr;

        if (mpidir == 2) {
            oxi = old_matriz->totalX;
            oyi = old_matriz->totalY;
            ozi = old_matriz->totalZ;

            this_obj.matriz.totalX = ozi;
            this_obj.matriz.totalY = oxi;
            this_obj.matriz.totalZ = oyi;

            oxi = old_despl->nX;
            oyi = old_despl->nY;
            ozi = old_despl->nZ;

            this_obj.despl.nX = ozi;
            this_obj.despl.nY = oxi;
            this_obj.despl.nZ = oyi;

            oxi = old_despl->mX1;
            oyi = old_despl->mY1;
            ozi = old_despl->mZ1;

            this_obj.despl.mX1 = ozi;
            this_obj.despl.mY1 = oxi;
            this_obj.despl.mZ1 = oyi;

            oxi = old_despl->mX2;
            oyi = old_despl->mY2;
            ozi = old_despl->mZ2;

            this_obj.despl.mX2 = ozi;
            this_obj.despl.mY2 = oxi;
            this_obj.despl.mZ2 = oyi;

            roxi = old_despl->originX;
            royi = old_despl->originY;
            rozi = old_despl->originZ;

            this_obj.despl.originX = rozi;
            this_obj.despl.originY = roxi;
            this_obj.despl.originZ = royi;

            poxi = &old_despl->desX;
            poyi = &old_despl->desY;
            pozi = &old_despl->desZ;

            this_obj.despl.desX = *pozi;
            this_obj.despl.desY = *poxi;
            this_obj.despl.desZ = *poyi;
        } else if (mpidir == 1) {
            oxi = old_matriz->totalX;
            oyi = old_matriz->totalY;
            ozi = old_matriz->totalZ;

            this_obj.matriz.totalX = oyi;
            this_obj.matriz.totalY = ozi;
            this_obj.matriz.totalZ = oxi;

            oxi = old_despl->nX;
            oyi = old_despl->nY;
            ozi = old_despl->nZ;

            this_obj.despl.nX = oyi;
            this_obj.despl.nY = ozi;
            this_obj.despl.nZ = oxi;

            oxi = old_despl->mX1;
            oyi = old_despl->mY1;
            ozi = old_despl->mZ1;

            this_obj.despl.mX1 = oyi;
            this_obj.despl.mY1 = ozi;
            this_obj.despl.mZ1 = oxi;

            oxi = old_despl->mX2;
            oyi = old_despl->mY2;
            ozi = old_despl->mZ2;

            this_obj.despl.mX2 = oyi;
            this_obj.despl.mY2 = ozi;
            this_obj.despl.mZ2 = oxi;

            roxi = old_despl->originX;
            royi = old_despl->originY;
            rozi = old_despl->originZ;

            this_obj.despl.originX = royi;
            this_obj.despl.originY = rozi;
            this_obj.despl.originZ = roxi;

            poxi = &old_despl->desX;
            poyi = &old_despl->desY;
            pozi = &old_despl->desZ;

            this_obj.despl.desX = *poyi;
            this_obj.despl.desY = *pozi;
            this_obj.despl.desZ = *poxi;
        }

        delete old_despl;
        delete old_matriz;
    }

    void rotate_generateCurrent_Field_Sources(Parseador_t& this_obj, int mpidir) {
        int tama = this_obj.nodsrc.n_nodSrc;
        for (int i = 0; i < tama; ++i) {
            int tama2 = this_obj.nodsrc.NodalSource[i].n_c1P;
            for (int ii = 0; ii < tama2; ++ii) {
                ROTATEMPI_SCALED(mpidir, this_obj.nodsrc.NodalSource[i].c1P[ii]);
            }
            int tama3 = this_obj.nodsrc.NodalSource[i].n_c2P;
            for (int ii = 0; ii < tama3; ++ii) {
                ROTATEMPI_SCALED(mpidir, this_obj.nodsrc.NodalSource[i].c2P[ii]);
            }
        }
    }

    void rotate_generatePlaneWaves(Parseador_t& this_obj, int mpidir) {
        PlaneWaves_t* old_plnSrc = new PlaneWaves_t(this_obj.plnSrc);
        int tama = this_obj.plnSrc.nc;
        for (int i = 0; i < tama; ++i) {
            double theta = old_plnSrc->collection[i].theta;
            double phi = old_plnSrc->collection[i].phi;
            double alpha = old_plnSrc->collection[i].alpha;
            double beta = old_plnSrc->collection[i].beta;

            if (mpidir == 2) {
                int oxi = old_plnSrc->collection[i].coor1[0];
                int oxe = old_plnSrc->collection[i].coor2[0];
                int oyi = old_plnSrc->collection[i].coor1[1];
                int oye = old_plnSrc->collection[i].coor2[1];
                int ozi = old_plnSrc->collection[i].coor1[2];
                int oze = old_plnSrc->collection[i].coor2[2];

                this_obj.plnSrc.collection[i].coor1[0] = ozi;
                this_obj.plnSrc.collection[i].coor2[0] = oze;
                this_obj.plnSrc.collection[i].coor1[1] = oxi;
                this_obj.plnSrc.collection[i].coor2[1] = oxe;
                this_obj.plnSrc.collection[i].coor1[2] = oyi;
                this_obj.plnSrc.collection[i].coor2[2] = oye;

                this_obj.plnSrc.collection[i].theta = std::atan2(std::sqrt(std::cos(theta) * std::cos(theta) + std::cos(phi) * std::cos(phi) * std::sin(theta) * std::sin(theta)), std::sin(phi) * std::sin(theta));
                this_obj.plnSrc.collection[i].phi = std::atan2(std::cos(phi) * std::sin(theta), std::cos(theta));
                this_obj.plnSrc.collection[i].alpha = std::atan2(std::sqrt(std::cos(alpha) * std::cos(alpha) + std::cos(beta) * std::cos(beta) * std::sin(alpha) * std::sin(alpha)), std::sin(beta) * std::sin(alpha));
                this_obj.plnSrc.collection[i].beta = std::atan2(std::cos(beta) * std::sin(alpha), std::cos(alpha));
            } else if (mpidir == 1) {
                int oxi = old_plnSrc->collection[i].coor1[0];
                int oxe = old_plnSrc->collection[i].coor2[0];
                int oyi = old_plnSrc->collection[i].coor1[1];
                int oye = old_plnSrc->collection[i].coor2[1];
                int ozi = old_plnSrc->collection[i].coor1[2];
                int oze = old_plnSrc->collection[i].coor2[2];

                this_obj.plnSrc.collection[i].coor1[0] = oyi;
                this_obj.plnSrc.collection[i].coor2[0] = oye;
                this_obj.plnSrc.collection[i].coor1[1] = ozi;
                this_obj.plnSrc.collection[i].coor2[1] = oze;
                this_obj.plnSrc.collection[i].coor1[2] = oxi;
                this_obj.plnSrc.collection[i].coor2[2] = oxe;

                this_obj.plnSrc.collection[i].theta = std::atan2(std::sqrt(std::cos(theta) * std::cos(theta) + std::sin(phi) * std::sin(phi) * std::sin(theta) * std::sin(theta)), std::cos(phi) * std::sin(theta));
                this_obj.plnSrc.collection[i].phi = std::atan2(std::cos(theta), std::sin(phi) * std::sin(theta));
                this_obj.plnSrc.collection[i].alpha = std::atan2(std::sqrt(std::cos(alpha) * std::cos(alpha) + std::sin(beta) * std::sin(beta) * std::sin(alpha) * std::sin(alpha)), std::cos(beta) * std::sin(alpha));
                this_obj.plnSrc.collection[i].beta = std::atan2(std::cos(alpha), std::sin(beta) * std::sin(alpha));
            }
        }
        delete old_plnSrc;
    }

    void rotate_generateBoxSources(Parseador_t& this_obj, int mpidir) {
        Boxes_t* old_boxSrc = new Boxes_t(this_obj.boxSrc);
        int tama = this_obj.boxSrc.nvols;
        for (int i = 0; i < tama; ++i) {
            if (mpidir == 2) {
                int oxi = old_boxSrc->vols[i].coor1[0];
                int oxe = old_boxSrc->vols[i].coor2[0];
                int oyi = old_boxSrc->vols[i].coor1[1];
                int oye = old_boxSrc->vols[i].coor2[1];
                int ozi = old_boxSrc->vols[i].coor1[2];
                int oze = old_boxSrc->vols[i].coor2[2];

                this_obj.boxSrc.vols[i].coor1[0] = ozi;
                this_obj.boxSrc.vols[i].coor2[0] = oze;
                this_obj.boxSrc.vols[i].coor1[1] = oxi;
                this_obj.boxSrc.vols[i].coor2[1] = oxe;
                this_obj.boxSrc.vols[i].coor1[2] = oyi;
                this_obj.boxSrc.vols[i].coor2[2] = oye;
            } else if (mpidir == 1) {
                int oxi = old_boxSrc->vols[i].coor1[0];
                int oxe = old_boxSrc->vols[i].coor2[0];
                int oyi = old_boxSrc->vols[i].coor1[1];
                int oye = old_boxSrc->vols[i].coor2[1];
                int ozi = old_boxSrc->vols[i].coor1[2];
                int oze = old_boxSrc->vols[i].coor2[2];

                this_obj.boxSrc.vols[i].coor1[0] = oyi;
                this_obj.boxSrc.vols[i].coor2[0] = oye;
                this_obj.boxSrc.vols[i].coor1[1] = ozi;
                this_obj.boxSrc.vols[i].coor2[1] = oze;
                this_obj.boxSrc.vols[i].coor1[2] = oxi;
                this_obj.boxSrc.vols[i].coor2[2] = oxe;
            }
        }
        delete old_boxSrc;
    }

    void rotate_generateFronteras(Parseador_t& this_obj, int mpidir) {
        int oxl, oxu, oyl, oyu, ozl, ozu;
        FronteraPML_t OPML_XL, OPML_XU, OPML_YL, OPML_YU, OPML_ZL, OPML_ZU;

        if (mpidir == 2) {
            oxl = this_obj.front.tipofrontera[0];
            oxu = this_obj.front.tipofrontera[1];
            oyl = this_obj.front.tipofrontera[2];
            oyu = this_obj.front.tipofrontera[3];
            ozl = this_obj.front.tipofrontera[4];
            ozu = this_obj.front.tipofrontera[5];

            this_obj.front.tipofrontera[0] = ozl;
            this_obj.front.tipofrontera[1] = ozu;
            this_obj.front.tipofrontera[2] = oxl;
            this_obj.front.tipofrontera[3] = oxu;
            this_obj.front.tipofrontera[4] = oyl;
            this_obj.front.tipofrontera[5] = oyu;

            OPML_XL.orden = this_obj.front.propiedadesPML[0].orden;
            OPML_XU.orden = this_obj.front.propiedadesPML[1].orden;
            OPML_YL.orden = this_obj.front.propiedadesPML[2].orden;
            OPML_YU.orden = this_obj.front.propiedadesPML[3].orden;
            OPML_ZL.orden = this_obj.front.propiedadesPML[4].orden;
            OPML_ZU.orden = this_obj.front.propiedadesPML[5].orden;

            this_obj.front.propiedadesPML[0].orden = OPML_ZL.orden;
            this_obj.front.propiedadesPML[1].orden = OPML_ZU.orden;
            this_obj.front.propiedadesPML[2].orden = OPML_XL.orden;
            this_obj.front.propiedadesPML[3].orden = OPML_XU.orden;
            this_obj.front.propiedadesPML[4].orden = OPML_YL.orden;
            this_obj.front.propiedadesPML[5].orden = OPML_YU.orden;

            OPML_XL.refl = this_obj.front.propiedadesPML[0].refl;
            OPML_XU.refl = this_obj.front.propiedadesPML[1].refl;
            OPML_YL.refl = this_obj.front.propiedadesPML[2].refl;
            OPML_YU.refl = this_obj.front.propiedadesPML[3].refl;
            OPML_ZL.refl = this_obj.front.propiedadesPML[4].refl;
            OPML_ZU.refl = this_obj.front.propiedadesPML[5].refl;

            this_obj.front.propiedadesPML[0].refl = OPML_ZL.refl;
            this_obj.front.propiedadesPML[1].refl = OPML_ZU.refl;
            this_obj.front.propiedadesPML[2].refl = OPML_XL.refl;
            this_obj.front.propiedadesPML[3].refl = OPML_XU.refl;
            this_obj.front.propiedadesPML[4].refl = OPML_YL.refl;
            this_obj.front.propiedadesPML[5].refl = OPML_YU.refl;

            OPML_XL.numCapas = this_obj.front.propiedadesPML[0].numCapas;
            OPML_XU.numCapas = this_obj.front.propiedadesPML[1].numCapas;
            OPML_YL.numCapas = this_obj.front.propiedadesPML[2].numCapas;
            OPML_YU.numCapas = this_obj.front.propiedadesPML[3].numCapas;
            OPML_ZL.numCapas = this_obj.front.propiedadesPML[4].numCapas;
            OPML_ZU.numCapas = this_obj.front.propiedadesPML[5].numCapas;

            this_obj.front.propiedadesPML[0].numCapas = OPML_ZL.numCapas;
            this_obj.front.propiedadesPML[1].numCapas = OPML_ZU.numCapas;
            this_obj.front.propiedadesPML[2].numCapas = OPML_XL.numCapas;
            this_obj.front.propiedadesPML[3].numCapas = OPML_XU.numCapas;
            this_obj.front.propiedadesPML[4].numCapas = OPML_YL.numCapas;
            this_obj.front.propiedadesPML[5].numCapas = OPML_YU.numCapas;
        } else if (mpidir == 1) {
            oxl = this_obj.front.tipofrontera[0];
            oxu = this_obj.front.tipofrontera[1];
            oyl = this_obj.front.tipofrontera[2];
            oyu = this_obj.front.tipofrontera[3];
            ozl = this_obj.front.tipofrontera[4];
            ozu = this_obj.front.tipofrontera[5];

            this_obj.front.tipofrontera[0] = oyl;
            this_obj.front.tipofrontera[1] = oyu;
            this_obj.front.tipofrontera[2] = ozl;
            this_obj.front.tipofrontera[3] = ozu;
            this_obj.front.tipofrontera[4] = oxl;
            this_obj.front.tipofrontera[5] = oxu;

            OPML_XL.orden = this_obj.front.propiedadesPML[0].orden;
            OPML_XU.orden = this_obj.front.propiedadesPML[1].orden;
            OPML_YL.orden = this_obj.front.propiedadesPML[2].orden;
            OPML_YU.orden = this_obj.front.propiedadesPML[3].orden;
            OPML_ZL.orden = this_obj.front.propiedadesPML[4].orden;
            OPML_ZU.orden = this_obj.front.propiedadesPML[5].orden;

            this_obj.front.propiedadesPML[0].orden = OPML_YL.orden;
            this_obj.front.propiedadesPML[1].orden = OPML_YU.orden;
            this_obj.front.propiedadesPML[2].orden = OPML_ZL.orden;
            this_obj.front.propiedadesPML[3].orden = OPML_ZU.orden;
            this_obj.front.propiedadesPML[4].orden = OPML_XL.orden;
            this_obj.front.propiedadesPML[5].orden = OPML_XU.orden;

            OPML_XL.refl = this_obj.front.propiedadesPML[0].refl;
            OPML_XU.refl = this_obj.front.propiedadesPML[1].refl;
            OPML_YL.refl = this_obj.front.propiedadesPML[2].refl;
            OPML_YU.refl = this_obj.front.propiedadesPML[3].refl;
            OPML_ZL.refl = this_obj.front.propiedadesPML[4].refl;
            OPML_ZU.refl = this_obj.front.propiedadesPML[5].refl;

            this_obj.front.propiedadesPML[0].refl = OPML_YL.refl;
            this_obj.front.propiedadesPML[1].refl = OPML_YU.refl;
            this_obj.front.propiedadesPML[2].refl = OPML_ZL.refl;
            this_obj.front.propiedadesPML[3].refl = OPML_ZU.refl;
            this_obj.front.propiedadesPML[4].refl = OPML_XL.refl;
            this_obj.front.propiedadesPML[5].refl = OPML_XU.refl;

            OPML_XL.numCapas = this_obj.front.propiedadesPML[0].numCapas;
            OPML_XU.numCapas = this_obj.front.propiedadesPML[1].numCapas;
            OPML_YL.numCapas = this_obj.front.propiedadesPML[2].numCapas;
            OPML_YU.numCapas = this_obj.front.propiedadesPML[3].numCapas;
            OPML_ZL.numCapas = this_obj.front.propiedadesPML[4].numCapas;
            OPML_ZU.numCapas = this_obj.front.propiedadesPML[5].numCapas;

            this_obj.front.propiedadesPML[0].numCapas = OPML_YL.numCapas;
            this_obj.front.propiedadesPML[1].numCapas = OPML_YU.numCapas;
            this_obj.front.propiedadesPML[2].numCapas = OPML_ZL.numCapas;
            this_obj.front.propiedadesPML[3].numCapas = OPML_ZU.numCapas;
            this_obj.front.propiedadesPML[4].numCapas = OPML_XL.numCapas;
            this_obj.front.propiedadesPML[5].numCapas = OPML_XU.numCapas;
        }
    }

    void rotate_generatePECs(Parseador_t& this_obj, int mpidir) {
        int tama = this_obj.pecregs.nvols;
        for (int i = 0; i < tama; ++i) {
            ROTATEMPI(mpidir, this_obj.pecRegs.Vols[i]);
        }
        tama = this_obj.pecregs.nsurfs;
        for (int i = 0; i < tama; ++i) {
            ROTATEMPI(mpidir, this_obj.pecRegs.Surfs[i]);
        }
        tama = this_obj.pecregs.nlins;
        for (int i = 0; i < tama; ++i) {
            ROTATEMPI(mpidir, this_obj.pecRegs.Lins[i]);
        }
    }

    void rotate_generatePMCs(Parseador_t& this_obj, int mpidir) {
        int tama = this_obj.pmcregs.nvols;
        for (int i = 0; i < tama; ++i) {
            ROTATEMPI(mpidir, this_obj.pmcRegs.Vols[i]);
        }
        tama = this_obj.pmcregs.nsurfs;
        for (int i = 0; i < tama; ++i) {
            ROTATEMPI(mpidir, this_obj.pmcRegs.Surfs[i]);
        }
        tama = this_obj.pmcregs.nlins;
        for (int i = 0; i < tama; ++i) {
            ROTATEMPI(mpidir, this_obj.pmcRegs.Lins[i]);
        }
    }

    void rotate_generateNONMetals(Parseador_t& this_obj, int mpidir) {
        int tama, tama2, tama3, i, ii;

        tama = this_obj.DielRegs.nvols;
        for (i = 0; i < tama; ++i) {
            tama2 = this_obj.DielRegs.vols[i].n_c1P;
            for (ii = 0; ii < tama2; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.vols[i].C1P[ii]);
            }
            if (tama2 > 0) {
                this_obj.DielRegs.vols[i].DiodOrI = this_obj.DielRegs.vols[i].c1P[tama2 - 1].Or;
            }
            tama3 = this_obj.DielRegs.vols[i].n_c2P;
            for (ii = 0; ii < tama3; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.vols[i].C2P[ii]);
            }
            if (tama3 > 0) {
                this_obj.DielRegs.vols[i].DiodOrI = this_obj.DielRegs.vols[i].c2P[tama3 - 1].Or;
            }
        }

        tama = this_obj.DielRegs.nsurfs;
        for (i = 0; i < tama; ++i) {
            tama2 = this_obj.DielRegs.surfs[i].n_c1P;
            for (ii = 0; ii < tama2; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.surfs[i].C1P[ii]);
            }
            if (tama2 > 0) {
                this_obj.DielRegs.surfs[i].DiodOrI = this_obj.DielRegs.surfs[i].c1P[tama2 - 1].Or;
            }
            tama3 = this_obj.DielRegs.surfs[i].n_c2P;
            for (ii = 0; ii < tama3; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.surfs[i].C2P[ii]);
            }
            if (tama3 > 0) {
                this_obj.DielRegs.surfs[i].DiodOrI = this_obj.DielRegs.surfs[i].c2P[tama3 - 1].Or;
            }
        }

        tama = this_obj.DielRegs.nlins;
        for (i = 0; i < tama; ++i) {
            tama2 = this_obj.DielRegs.lins[i].n_c1P;
            for (ii = 0; ii < tama2; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.lins[i].C1P[ii]);
            }
            if (tama2 > 0) {
                this_obj.DielRegs.lins[i].DiodOrI = this_obj.DielRegs.lins[i].c1P[tama2 - 1].Or;
            }
            tama3 = this_obj.DielRegs.lins[i].n_c2P;
            for (ii = 0; ii < tama3; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.lins[i].C2P[ii]);
            }
            if (tama3 > 0) {
                this_obj.DielRegs.lins[i].DiodOrI = this_obj.DielRegs.lins[i].c2P[tama3 - 1].Or;
            }
        }
    }

    void rotate_generateANISOTROPICs(Parseador_t& this_obj, int mpidir) {
        if ((mpidir != 1) && (this_obj.ANIMATS.nvols + this_obj.ANIMATS.nsurfs + this_obj.ANIMATS.nlins != 0)) {
            std::cout << "Rotations in anisotropic unsupported" << std::endl;
            std::exit(1);
        }
    }

    void rotate_generateThinWires(Parseador_t& this_obj, int mpidir) {
        int tama = this_obj.twires.n_tw;
        for (int i = 0; i < tama; ++i) {
            int tama2 = this_obj.twires.TW[i].N_TWC;
            for (int ii = 0; ii < tama2; ++ii) {
                int oldx = this_obj.twires.tw[i].tWc[ii].i;
                int oldy = this_obj.twires.tw[i].tWc[ii].j;
                int oldz = this_obj.twires.tw[i].tWc[ii].K;

                if (mpidir == 2) {
                    this_obj.twires.tw[i].tWc[ii].i = oldz;
                    this_obj.twires.tw[i].tWc[ii].j = oldx;
                    this_obj.twires.tw[i].tWc[ii].K = oldy;

                    switch (this_obj.twires.tw[i].tWc[ii].d) {
                        case iEx:
                            this_obj.twires.tw[i].tWc[ii].d = iEy;
                            break;
                        case iEy:
                            this_obj.twires.tw[i].tWc[ii].d = iEz;
                            break;
                        case iEz:
                            this_obj.twires.tw[i].tWc[ii].d = iEx;
                            break;
                    }
                } else if (mpidir == 1) {
                    this_obj.twires.tw[i].tWc[ii].i = oldy;
                    this_obj.twires.tw[i].tWc[ii].j = oldz;
                    this_obj.twires.tw[i].tWc[ii].K = oldx;

                    switch (this_obj.twires.tw[i].tWc[ii].d) {
                        case iEx:
                            this_obj.twires.tw[i].tWc[ii].d = i