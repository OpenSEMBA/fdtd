#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>

// Forward declarations and includes for types defined in NFDETypes_m
// Assuming NFDETypes_m defines the following structs/classes. 
// Since the full definition isn't provided, we assume standard mappings based on usage.

// Placeholder for types from NFDETypes_m to ensure compilation context
// In a real scenario, these would be included from a header file.

struct Parseador_t;
struct Desplazamiento_t;
struct MatrizMedios_t;
struct mtln_t;
struct PlaneWaves_t;
struct Boxes_t;
struct FronteraPML_t;
struct NodalSource_t;
struct PEC_Regs_t;
struct PMC_Regs_t;
struct DielRegs_t;
struct ThinWires_t;
struct SlantedWiresInfo_t;

// Constants and Types assumed from NFDETypes_m
using RK = double;
using RKIND = double;

// Enumerations assumed from context
enum Direction { iEx, iEy, iEz };

// Helper function for sign
inline int sign(int val, int base) {
    return (base < 0) ? -val : val;
}

// Helper function for MPI rotation of coordinates (assumed implementation based on usage)
// Note: The actual implementation of ROTATEMPI and ROTATEMPI_SCALED is not provided in the snippet.
// We provide a stub that performs the coordinate permutation logic seen in other functions.
void ROTATEMPI(int mpidir, std::vector<double>& coords) {
    if (coords.size() < 3) return;
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    
    if (mpidir == 2) {
        // X->Y->Z->X
        coords[0] = z;
        coords[1] = x;
        coords[2] = y;
    } else if (mpidir == 1) {
        // X->Z->Y->X
        coords[0] = y;
        coords[1] = z;
        coords[2] = x;
    }
}

void ROTATEMPI_SCALED(int mpidir, std::vector<double>& coords) {
    ROTATEMPI(mpidir, coords);
}

namespace nfde_rotate_m {

    void nfde_rotate(Parseador_t& this_obj, int mpidir) {
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
        // Deep copy of mtln
        mtln_t old_mtln = this_obj.mtln; 
        
        int tama_cables = old_mtln.cables.size();
        for (int i = 0; i < tama_cables; ++i) {
            // Assuming cables(i)%ptr%segments is a vector of segments
            // If ptr is a shared_ptr or raw pointer, we access the object it points to.
            // For simplicity, assuming direct access or that 'cables' contains the data directly in C++ translation
            // If 'ptr' implies indirection, we dereference.
            // Let's assume old_mtln.cables[i] has a 'segments' member directly for the C++ struct mapping
            // or old_mtln.cables[i].ptr->segments.
            
            // To be safe with the "pointer" semantics in Fortran, we assume the C++ struct mirrors the data layout
            // but uses vectors.
            
            auto& segments = old_mtln.cables[i].segments; 
            int tama_segments = segments.size();
            for (int j = 0; j < tama_segments; ++j) {
                int x = segments[j].x;
                int y = segments[j].y;
                int z = segments[j].z;
                int or_val = segments[j].orientation;
                
                if (mpidir == 2) {
                    this_obj.mtln.cables[i].segments[j].x = z;
                    this_obj.mtln.cables[i].segments[j].y = x;
                    this_obj.mtln.cables[i].segments[j].z = y;
                    
                    int abs_or = std::abs(or_val);
                    if (abs_or == 1) {
                        this_obj.mtln.cables[i].segments[j].orientation = sign(2, or_val);
                    } else if (abs_or == 2) {
                        this_obj.mtln.cables[i].segments[j].orientation = sign(3, or_val);
                    } else if (abs_or == 3) {
                        this_obj.mtln.cables[i].segments[j].orientation = sign(1, or_val);
                    }
                } else if (mpidir == 1) {
                    this_obj.mtln.cables[i].segments[j].x = y;
                    this_obj.mtln.cables[i].segments[j].y = z;
                    this_obj.mtln.cables[i].segments[j].z = x;
                    
                    int abs_or = std::abs(or_val);
                    if (abs_or == 1) {
                        this_obj.mtln.cables[i].segments[j].orientation = sign(3, or_val);
                    } else if (abs_or == 2) {
                        this_obj.mtln.cables[i].segments[j].orientation = sign(1, or_val);
                    } else if (abs_or == 3) {
                        this_obj.mtln.cables[i].segments[j].orientation = sign(2, or_val);
                    }
                }
            }
        }
    }
#endif

    void rotate_generateSpaceSteps(Parseador_t& this_obj, int mpidir) {
        Desplazamiento_t old_despl = this_obj.despl;
        MatrizMedios_t old_matriz = this_obj.matriz;

        int oxi, oyi, ozi;
        double roxi, royi, rozi;
        
        // Note: In Fortran, pointers poxi, poyi, pozi point to arrays inside old_despl.
        // In C++, we copy the vectors.
        std::vector<double> poxi = old_despl.desX;
        std::vector<double> poyi = old_despl.desY;
        std::vector<double> pozi = old_despl.desZ;

        if (mpidir == 2) {
            // X->Y->Z->X
            oxi = old_matriz.totalX;
            oyi = old_matriz.totalY;
            ozi = old_matriz.totalZ;
            
            this_obj.matriz.totalX = ozi;
            this_obj.matriz.totalY = oxi;
            this_obj.matriz.totalZ = oyi;
            
            oxi = old_despl.nX;
            oyi = old_despl.nY;
            ozi = old_despl.nZ;
            
            this_obj.despl.nX = ozi;
            this_obj.despl.nY = oxi;
            this_obj.despl.nZ = oyi;
            
            oxi = old_despl.mX1;
            oyi = old_despl.mY1;
            ozi = old_despl.mZ1;
            
            this_obj.despl.mX1 = ozi;
            this_obj.despl.mY1 = oxi;
            this_obj.despl.mZ1 = oyi;
            
            oxi = old_despl.mX2;
            oyi = old_despl.mY2;
            ozi = old_despl.mZ2;
            
            this_obj.despl.mX2 = ozi;
            this_obj.despl.mY2 = oxi;
            this_obj.despl.mZ2 = oyi;
            
            roxi = old_despl.originX;
            royi = old_despl.originY;
            rozi = old_despl.originZ;
            
            this_obj.despl.originX = rozi;
            this_obj.despl.originY = roxi;
            this_obj.despl.originZ = royi;
            
            // Pointer assignment in Fortran becomes vector assignment in C++
            this_obj.despl.desX = pozi;
            this_obj.despl.desY = poxi;
            this_obj.despl.desZ = poyi;
            
        } else if (mpidir == 1) {
            // X->Z->Y->X
            oxi = old_matriz.totalX;
            oyi = old_matriz.totalY;
            ozi = old_matriz.totalZ;
            
            this_obj.matriz.totalX = oyi;
            this_obj.matriz.totalY = ozi;
            this_obj.matriz.totalZ = oxi;
            
            oxi = old_despl.nX;
            oyi = old_despl.nY;
            ozi = old_despl.nZ;
            
            this_obj.despl.nX = oyi;
            this_obj.despl.nY = ozi;
            this_obj.despl.nZ = oxi;
            
            oxi = old_despl.mX1;
            oyi = old_despl.mY1;
            ozi = old_despl.mZ1;
            
            this_obj.despl.mX1 = oyi;
            this_obj.despl.mY1 = ozi;
            this_obj.despl.mZ1 = oxi;
            
            oxi = old_despl.mX2;
            oyi = old_despl.mY2;
            ozi = old_despl.mZ2;
            
            this_obj.despl.mX2 = oyi;
            this_obj.despl.mY2 = ozi;
            this_obj.despl.mZ2 = oxi;
            
            roxi = old_despl.originX;
            royi = old_despl.originY;
            rozi = old_despl.originZ;
            
            this_obj.despl.originX = royi;
            this_obj.despl.originY = rozi;
            this_obj.despl.originZ = roxi;
            
            this_obj.despl.desX = poyi;
            this_obj.despl.desY = pozi;
            this_obj.despl.desZ = poxi;
        }
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
        PlaneWaves_t old_plnSrc = this_obj.plnSrc;
        int tama = this_obj.plnSrc.nc;
        
        for (int i = 0; i < tama; ++i) {
            double theta = old_plnSrc.collection[i].theta;
            double phi = old_plnSrc.collection[i].phi;
            double alpha = old_plnSrc.collection[i].alpha;
            double beta = old_plnSrc.collection[i].beta;
            
            int oxi, oxe, oyi, oye, ozi, oze;
            
            if (mpidir == 2) {
                oxi = old_plnSrc.collection[i].coor1[0];
                oxe = old_plnSrc.collection[i].coor2[0];
                oyi = old_plnSrc.collection[i].coor1[1];
                oye = old_plnSrc.collection[i].coor2[1];
                ozi = old_plnSrc.collection[i].coor1[2];
                oze = old_plnSrc.collection[i].coor2[2];
                
                this_obj.plnSrc.collection[i].coor1[0] = ozi;
                this_obj.plnSrc.collection[i].coor2[0] = oze;
                this_obj.plnSrc.collection[i].coor1[1] = oxi;
                this_obj.plnSrc.collection[i].coor2[1] = oxe;
                this_obj.plnSrc.collection[i].coor1[2] = oyi;
                this_obj.plnSrc.collection[i].coor2[2] = oye;
                
                this_obj.plnSrc.collection[i].theta = std::atan2(std::sqrt(std::cos(theta)*std::cos(theta) + std::cos(phi)*std::cos(phi)*std::sin(theta)*std::sin(theta)), std::sin(phi)*std::sin(theta));
                this_obj.plnSrc.collection[i].phi = std::atan2(std::cos(phi)*std::sin(theta), std::cos(theta));
                this_obj.plnSrc.collection[i].alpha = std::atan2(std::sqrt(std::cos(alpha)*std::cos(alpha) + std::cos(beta)*std::cos(beta)*std::sin(alpha)*std::sin(alpha)), std::sin(beta)*std::sin(alpha));
                this_obj.plnSrc.collection[i].beta = std::atan2(std::cos(beta)*std::sin(alpha), std::cos(alpha));
                
            } else if (mpidir == 1) {
                oxi = old_plnSrc.collection[i].coor1[0];
                oxe = old_plnSrc.collection[i].coor2[0];
                oyi = old_plnSrc.collection[i].coor1[1];
                oye = old_plnSrc.collection[i].coor2[1];
                ozi = old_plnSrc.collection[i].coor1[2];
                oze = old_plnSrc.collection[i].coor2[2];
                
                this_obj.plnSrc.collection[i].coor1[0] = oyi;
                this_obj.plnSrc.collection[i].coor2[0] = oye;
                this_obj.plnSrc.collection[i].coor1[1] = ozi;
                this_obj.plnSrc.collection[i].coor2[1] = oze;
                this_obj.plnSrc.collection[i].coor1[2] = oxi;
                this_obj.plnSrc.collection[i].coor2[2] = oxe;
                
                this_obj.plnSrc.collection[i].theta = std::atan2(std::sqrt(std::cos(theta)*std::cos(theta) + std::sin(phi)*std::sin(phi)*std::sin(theta)*std::sin(theta)), std::cos(phi)*std::sin(theta));
                this_obj.plnSrc.collection[i].phi = std::atan2(std::cos(theta), std::sin(phi)*std::sin(theta));
                this_obj.plnSrc.collection[i].alpha = std::atan2(std::sqrt(std::cos(alpha)*std::cos(alpha) + std::sin(beta)*std::sin(beta)*std::sin(alpha)*std::sin(alpha)), std::cos(beta)*std::sin(alpha));
                this_obj.plnSrc.collection[i].beta = std::atan2(std::cos(alpha), std::sin(beta)*std::sin(alpha));
            }
        }
    }

    void rotate_generateBoxSources(Parseador_t& this_obj, int mpidir) {
        Boxes_t old_boxSrc = this_obj.boxSrc;
        int tama = this_obj.boxSrc.nvols;
        
        for (int i = 0; i < tama; ++i) {
            int oxi, oxe, oyi, oye, ozi, oze;
            
            if (mpidir == 2) {
                oxi = old_boxSrc.vols[i].coor1[0];
                oxe = old_boxSrc.vols[i].coor2[0];
                oyi = old_boxSrc.vols[i].coor1[1];
                oye = old_boxSrc.vols[i].coor2[1];
                ozi = old_boxSrc.vols[i].coor1[2];
                oze = old_boxSrc.vols[i].coor2[2];
                
                this_obj.boxSrc.vols[i].coor1[0] = ozi;
                this_obj.boxSrc.vols[i].coor2[0] = oze;
                this_obj.boxSrc.vols[i].coor1[1] = oxi;
                this_obj.boxSrc.vols[i].coor2[1] = oxe;
                this_obj.boxSrc.vols[i].coor1[2] = oyi;
                this_obj.boxSrc.vols[i].coor2[2] = oye;
                
            } else if (mpidir == 1) {
                oxi = old_boxSrc.vols[i].coor1[0];
                oxe = old_boxSrc.vols[i].coor2[0];
                oyi = old_boxSrc.vols[i].coor1[1];
                oye = old_boxSrc.vols[i].coor2[1];
                ozi = old_boxSrc.vols[i].coor1[2];
                oze = old_boxSrc.vols[i].coor2[2];
                
                this_obj.boxSrc.vols[i].coor1[0] = oyi;
                this_obj.boxSrc.vols[i].coor2[0] = oye;
                this_obj.boxSrc.vols[i].coor1[1] = ozi;
                this_obj.boxSrc.vols[i].coor2[1] = oze;
                this_obj.boxSrc.vols[i].coor1[2] = oxi;
                this_obj.boxSrc.vols[i].coor2[2] = oxe;
            }
        }
    }

    void rotate_generateFronteras(Parseador_t& this_obj, int mpidir) {
        int oxl, oxu, oyl, oyu, ozl, ozu;
        
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
            
            // PML Properties Rotation
            // Indices 0,1 are X, 2,3 are Y, 4,5 are Z in original array
            // After rotation X->Y->Z->X:
            // New X (0,1) comes from Old Z (4,5)
            // New Y (2,3) comes from Old X (0,1)
            // New Z (4,5) comes from Old Y (2,3)
            
            // Temp storage
            int orden_xl = this_obj.front.propiedadesPML[0].orden;
            int orden_xu = this_obj.front.propiedadesPML[1].orden;
            int orden_yl = this_obj.front.propiedadesPML[2].orden;
            int orden_yu = this_obj.front.propiedadesPML[3].orden;
            int orden_zl = this_obj.front.propiedadesPML[4].orden;
            int orden_zu = this_obj.front.propiedadesPML[5].orden;
            
            this_obj.front.propiedadesPML[0].orden = orden_zl;
            this_obj.front.propiedadesPML[1].orden = orden_zu;
            this_obj.front.propiedadesPML[2].orden = orden_xl;
            this_obj.front.propiedadesPML[3].orden = orden_xu;
            this_obj.front.propiedadesPML[4].orden = orden_yl;
            this_obj.front.propiedadesPML[5].orden = orden_yu;
            
            int refl_xl = this_obj.front.propiedadesPML[0].refl;
            int refl_xu = this_obj.front.propiedadesPML[1].refl;
            int refl_yl = this_obj.front.propiedadesPML[2].refl;
            int refl_yu = this_obj.front.propiedadesPML[3].refl;
            int refl_zl = this_obj.front.propiedadesPML[4].refl;
            int refl_zu = this_obj.front.propiedadesPML[5].refl;
            
            this_obj.front.propiedadesPML[0].refl = refl_zl;
            this_obj.front.propiedadesPML[1].refl = refl_zu;
            this_obj.front.propiedadesPML[2].refl = refl_xl;
            this_obj.front.propiedadesPML[3].refl = refl_xu;
            this_obj.front.propiedadesPML[4].refl = refl_yl;
            this_obj.front.propiedadesPML[5].refl = refl_yu;
            
            int numcapas_xl = this_obj.front.propiedadesPML[0].numCapas;
            int numcapas_xu = this_obj.front.propiedadesPML[1].numCapas;
            int numcapas_yl = this_obj.front.propiedadesPML[2].numCapas;
            int numcapas_yu = this_obj.front.propiedadesPML[3].numCapas;
            int numcapas_zl = this_obj.front.propiedadesPML[4].numCapas;
            int numcapas_zu = this_obj.front.propiedadesPML[5].numCapas;
            
            this_obj.front.propiedadesPML[0].numCapas = numcapas_zl;
            this_obj.front.propiedadesPML[1].numCapas = numcapas_zu;
            this_obj.front.propiedadesPML[2].numCapas = numcapas_xl;
            this_obj.front.propiedadesPML[3].numCapas = numcapas_xu;
            this_obj.front.propiedadesPML[4].numCapas = numcapas_yl;
            this_obj.front.propiedadesPML[5].numCapas = numcapas_yu;
            
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
            
            // X->Z->Y->X
            // New X (0,1) comes from Old Y (2,3)
            // New Y (2,3) comes from Old Z (4,5)
            // New Z (4,5) comes from Old X (0,1)
            
            int orden_xl = this_obj.front.propiedadesPML[0].orden;
            int orden_xu = this_obj.front.propiedadesPML[1].orden;
            int orden_yl = this_obj.front.propiedadesPML[2].orden;
            int orden_yu = this_obj.front.propiedadesPML[3].orden;
            int orden_zl = this_obj.front.propiedadesPML[4].orden;
            int orden_zu = this_obj.front.propiedadesPML[5].orden;
            
            this_obj.front.propiedadesPML[0].orden = orden_yl;
            this_obj.front.propiedadesPML[1].orden = orden_yu;
            this_obj.front.propiedadesPML[2].orden = orden_zl;
            this_obj.front.propiedadesPML[3].orden = orden_zu;
            this_obj.front.propiedadesPML[4].orden = orden_xl;
            this_obj.front.propiedadesPML[5].orden = orden_xu;
            
            int refl_xl = this_obj.front.propiedadesPML[0].refl;
            int refl_xu = this_obj.front.propiedadesPML[1].refl;
            int refl_yl = this_obj.front.propiedadesPML[2].refl;
            int refl_yu = this_obj.front.propiedadesPML[3].refl;
            int refl_zl = this_obj.front.propiedadesPML[4].refl;
            int refl_zu = this_obj.front.propiedadesPML[5].refl;
            
            this_obj.front.propiedadesPML[0].refl = refl_yl;
            this_obj.front.propiedadesPML[1].refl = refl_yu;
            this_obj.front.propiedadesPML[2].refl = refl_zl;
            this_obj.front.propiedadesPML[3].refl = refl_zu;
            this_obj.front.propiedadesPML[4].refl = refl_xl;
            this_obj.front.propiedadesPML[5].refl = refl_xu;
            
            int numcapas_xl = this_obj.front.propiedadesPML[0].numCapas;
            int numcapas_xu = this_obj.front.propiedadesPML[1].numCapas;
            int numcapas_yl = this_obj.front.propiedadesPML[2].numCapas;
            int numcapas_yu = this_obj.front.propiedadesPML[3].numCapas;
            int numcapas_zl = this_obj.front.propiedadesPML[4].numCapas;
            int numcapas_zu = this_obj.front.propiedadesPML[5].numCapas;
            
            this_obj.front.propiedadesPML[0].numCapas = numcapas_yl;
            this_obj.front.propiedadesPML[1].numCapas = numcapas_yu;
            this_obj.front.propiedadesPML[2].numCapas = numcapas_zl;
            this_obj.front.propiedadesPML[3].numCapas = numcapas_zu;
            this_obj.front.propiedadesPML[4].numCapas = numcapas_xl;
            this_obj.front.propiedadesPML[5].numCapas = numcapas_xu;
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
        // Volumes
        int tama = this_obj.DielRegs.nvols;
        for (int i = 0; i < tama; ++i) {
            int tama2 = this_obj.DielRegs.vols[i].n_c1P;
            for (int ii = 0; ii < tama2; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.vols[i].C1P[ii]);
            }
            if (tama2 > 0) {
                this_obj.DielRegs.vols[i].DiodOrI = this_obj.DielRegs.vols[i].c1P[tama2-1].Or;
            }
            int tama3 = this_obj.DielRegs.vols[i].n_c2P;
            for (int ii = 0; ii < tama3; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.vols[i].C2P[ii]);
            }
            if (tama3 > 0) {
                this_obj.DielRegs.vols[i].DiodOrI = this_obj.DielRegs.vols[i].c2P[tama3-1].Or;
            }
        }
        
        // Surfaces
        tama = this_obj.DielRegs.nsurfs;
        for (int i = 0; i < tama; ++i) {
            int tama2 = this_obj.DielRegs.surfs[i].n_c1P;
            for (int ii = 0; ii < tama2; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.surfs[i].C1P[ii]);
            }
            if (tama2 > 0) {
                this_obj.DielRegs.surfs[i].DiodOrI = this_obj.DielRegs.surfs[i].c1P[tama2-1].Or;
            }
            int tama3 = this_obj.DielRegs.surfs[i].n_c2P;
            for (int ii = 0; ii < tama3; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.surfs[i].C2P[ii]);
            }
            if (tama3 > 0) {
                this_obj.DielRegs.surfs[i].DiodOrI = this_obj.DielRegs.surfs[i].c2P[tama3-1].Or;
            }
        }
        
        // Lines
        tama = this_obj.DielRegs.nlins;
        for (int i = 0; i < tama; ++i) {
            int tama2 = this_obj.DielRegs.lins[i].n_c1P;
            for (int ii = 0; ii < tama2; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.lins[i].C1P[ii]);
            }
            if (tama2 > 0) {
                this_obj.DielRegs.lins[i].DiodOrI = this_obj.DielRegs.lins[i].c1P[tama2-1].Or;
            }
            int tama3 = this_obj.DielRegs.lins[i].n_c2P;
            for (int ii = 0; ii < tama3; ++ii) {
                ROTATEMPI(mpidir, this_obj.DielRegs.lins[i].C2P[ii]);
            }
            if (tama3 > 0) {
                this_obj.DielRegs.lins[i].DiodOrI = this_obj.DielRegs.lins[i].c2P[tama3-1].Or;
            }
        }
    }

    void rotate_generateANISOTROPICs(Parseador_t& this_obj, int mpidir) {
        if ((mpidir != 1) && ((this_obj.ANIMATS.nvols + this_obj.ANIMATS.nsurfs + this_obj.ANIMATS.nlins) != 0)) {
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
                        default:
                            break;
                    }
                } else if (mpidir == 1) {
                    this_obj.twires.tw[i].tWc[ii].i = oldy;
                    this_obj.twires.tw[i].tWc[ii].j = oldz;
                    this_obj.twires.tw[i].tWc[ii].K = oldx;
                    
                    switch (this_obj.twires.tw[i].tWc[ii].d) {
                        case iEx:
                            this_obj.twires.tw[i].tWc[ii].d = iEz;
                            break;
                        case iEy:
                            this_obj.twires.tw[i].tWc[ii].d = iEx;
                            break;
                        case iEz:
                            this_obj.twires.tw[i].tWc[ii].d = iEy;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }

    void rotate_generateSlantedWires(Parseador_t& this_obj, int mpidir) {
        int tama = this_obj.swires.n_sw;
        SlantedWiresInfo_t old_swires = this_obj.swires;
        
        for (int i = 0; i < tama; ++i) {
            int tama2 = this_obj.swires.sW[i].N_SWC;
            for (int ii = 0; ii < tama2; ++ii) {
                double oldx, oldy, oldz;
                
                if (mpidir == 2) {
                    oldx = this_obj.swires.sw[i].swc[ii].x;
                    oldy = this_obj.swires.sw[i].swc[ii].y;
                    oldz = this_obj.swires.sw[i].swc[ii].z;
                    
                    this_obj.swires.sw[i].swc[ii].x = oldz;
                    this_obj.swires.sw[i].swc[ii].y = oldx;
                    this_obj.swires.sw[i].swc[ii].z = oldy;
                } else if (mpidir == 1) {
                    oldx = this_obj.swires.sw[i].swc[ii].x;
                    oldy = this_obj.swires.sw[i].swc[ii].y;
                    oldz = this_obj.swires.sw[i].swc[ii].z;
                    
                    this_obj.swires.sw[i].swc[ii].x = oldy;
                    this_obj.swires.sw[i].swc[ii].y = oldz;
                    this_obj.swires.sw[i].swc[ii].z = oldx;
                }
            }
        }
    }

    // Stubs for remaining subroutines not fully defined in the input snippet
    void rotate_generateThinSlots(Parseador_t& this_obj, int mpidir) {
        // Implementation not provided in snippet
    }

    void rotate_generateLossyThinSurface(Parseador_t& this_obj, int mpidir) {
        // Implementation not provided in snippet
    }

    void rotate_generateFDMs(Parseador_t& this_obj, int mpidir) {
        // Implementation not provided in snippet
    }

    void rotate_generateSONDAs(Parseador_t& this_obj, int mpidir) {
        // Implementation not provided in snippet
    }

    void rotate_generateMasSondas(Parseador_t& this_obj, int mpidir) {
        // Implementation not provided in snippet
    }

    void rotate_generateBloqueProbes(Parseador_t& this_obj, int mpidir) {
        // Implementation not provided in snippet
    }

    void rotate_generateVolumicProbes(Parseador_t& this_obj, int mpidir) {
        // Implementation not provided in snippet
    }

} // namespace nfde_rotate_m

oldz = this.swires.sw[i].swc[ii].z;

                      this.swires.sw[i].swc[ii].x = oldz;
                      this.swires.sw[i].swc[ii].y = oldx;
                      this.swires.sw[i].swc[ii].z = oldy;
             } else if (MPIDIR == 1) {
                      oldx = this.swires.sw[i].swc[ii].x;
                      oldy = this.swires.sw[i].swc[ii].y;
                      oldz = this.swires.sw[i].swc[ii].z;

                      this.swires.sw[i].swc[ii].x = oldy;
                      this.swires.sw[i].swc[ii].y = oldz;
                      this.swires.sw[i].swc[ii].z = oldx;
             }
         }
         // FIN
      }
      allocate(old_swires);

      return;
   }

   void rotate_generateThinSlots(Parseador_t& this, int mpidir) {
      ThinSlots_t* old_tSlots = nullptr;
      int tama, tama2, i, ii;
      int oldx, oldy, oldz;

      tama = this.tSlots.n_Tg;
      old_tSlots = new ThinSlots_t(this.tSlots);
      for (i = 1; i <= tama; i++) {
         tama2 = this.tSlots.Tg[i].N_Tgc;
         for (ii = 1; ii <= tama2; ii++) {
              // ROTATE THIN SLOT
              if (MPIDIR == 2) {
                       oldx = this.tSlots.Tg[i].TgC[ii].i;
                       oldy = this.tSlots.Tg[i].TgC[ii].j;
                       oldz = this.tSlots.Tg[i].TgC[ii].K;

                       this.tSlots.Tg[i].TgC[ii].i = oldz;
                       this.tSlots.Tg[i].TgC[ii].j = oldx;
                       this.tSlots.Tg[i].TgC[ii].K = oldy;
                       switch (old_tSlots->Tg[i].TgC[ii].dir) {
                        case iEx:
                          this.tSlots.Tg[i].TgC[ii].dir = iEy;
                          break;
                        case iEY:
                          this.tSlots.Tg[i].TgC[ii].dir = iEz;
                          break;
                        case iEZ:
                          this.tSlots.Tg[i].TgC[ii].dir = iEx;
                          break;
                        default:
                          break;
                       }
              } else if (MPIDIR == 1) {
                       oldx = this.tSlots.Tg[i].TgC[ii].i;
                       oldy = this.tSlots.Tg[i].TgC[ii].j;
                       oldz = this.tSlots.Tg[i].TgC[ii].K;

                       this.tSlots.Tg[i].TgC[ii].i = oldy;
                       this.tSlots.Tg[i].TgC[ii].j = oldz;
                       this.tSlots.Tg[i].TgC[ii].K = oldx;

                       switch (old_tSlots->Tg[i].TgC[ii].dir) {
                        case iEx:
                          this.tSlots.Tg[i].TgC[ii].dir = iEz;
                          break;
                        case iEY:
                          this.tSlots.Tg[i].TgC[ii].dir = iEx;
                          break;
                        case iEZ:
                          this.tSlots.Tg[i].TgC[ii].dir = iEy;
                          break;
                        default:
                          break;
                       }
              }
         }
      }
      delete old_tSlots;
      return;
   }

   void rotate_generateLossyThinSurface(Parseador_t& this, int mpidir) {
      int tama2, tama, i, ii;

      tama = this.LossyThinSurfs.length;
      for (i = 1; i <= tama; i++) {
         tama2 = this.LossyThinSurfs.cs[i].nc;
         for (ii = 1; ii <= tama2; ii++) {
            ROTATEMPI(mpidir, this.LossyThinSurfs.cs[i].C[ii]);
         }
      }

      return;
   }

   void rotate_generateFDMs(Parseador_t& this, int mpidir) {
      int tama, tama2, i, ii;

      tama = (this.FRQDEPMATS.nvols);
      for (i = 1; i <= tama; i++) {
         tama2 = this.FRQDEPMATS.vols[i].n_C;
         for (ii = 1; ii <= tama2; ii++) {
            ROTATEMPI(mpidir, this.FRQDEPMATS.Vols[i].c[ii]);
         }
         rotate_freq_depend_material_properties(mpidir, this.FRQDEPMATS.Vols[i]);
      }

      tama = (this.FRQDEPMATS.nsurfs);
      for (i = 1; i <= tama; i++) {
         tama2 = this.FRQDEPMATS.surfs[i].n_C;
         for (ii = 1; ii <= tama2; ii++) {
            ROTATEMPI(mpidir, this.FRQDEPMATS.Surfs[i].c[ii]);
         }
         rotate_freq_depend_material_properties(mpidir, this.FRQDEPMATS.Surfs[i]);
      }

      tama = (this.FRQDEPMATS.nlins);
      for (i = 1; i <= tama; i++) {
         tama2 = this.FRQDEPMATS.Lins[i].n_C;
         for (ii = 1; ii <= tama2; ii++) {
            ROTATEMPI(mpidir, this.FRQDEPMATS.Lins[i].c[ii]);
         }
         rotate_freq_depend_material_properties(mpidir, this.FRQDEPMATS.Lins[i]);
      }
      return;
   }

   void rotate_generateSONDAs(Parseador_t& this, int mpidir) {
      int tama, tama2, tama3, i, ii, iii;
      FarField_Sonda_t* old_FarField = nullptr;
      Electric_Sonda_t* old_Electric = nullptr;
      Magnetic_Sonda_t* old_Magnetic = nullptr;
      double THETASTART, THETASTOP, PHISTART, PHISTOP;
      int iox, ioy, ioz;

      tama = this.oldSONDA.n_probes;
      // tres posibilidades FarField, Electric, Magnetic
      for (i = 1; i <= tama; i++) {
         tama2 = (this.oldSONDA.probes[i].n_FarField);
         for (ii = 1; ii <= tama2; ii++) {
            old_FarField = new FarField_Sonda_t(this.oldSONDA.probes[i].FarField[ii]);
            THETASTART = old_FarField->probe.thetastart;
            THETASTOP = old_FarField->probe.thetastop;
            PHISTART = old_FarField->probe.phistart;
            PHISTOP = old_FarField->probe.phistop;
            // mpirotate angulos farfield .... las coordenadas se rotan luego
            if (MPIDIR == 2) {
                   this.oldSONDA.probes[i].FarField[ii].probe.thetastart = atan2(sqrt(pow(cos(THETASTART), 2.0) + pow(cos(PHISTART), 2) * pow(sin(THETASTART), 2)), sin(PHISTART) * sin(THETASTART));
                   this.oldSONDA.probes[i].FarField[ii].probe.phistart = atan2(cos(PHISTART) * sin(THETASTART), cos(THETASTART));
                   this.oldSONDA.probes[i].FarField[ii].probe.thetastop = atan2(sqrt(pow(cos(THETASTOP), 2.0) + pow(cos(PHISTOP), 2) * pow(sin(THETASTOP), 2)), sin(PHISTOP) * sin(THETASTOP));
                   this.oldSONDA.probes[i].FarField[ii].probe.phistop = atan2(cos(PHISTOP) * sin(THETASTOP), cos(THETASTOP));
            } else if (MPIDIR == 1) {
                   this.oldSONDA.probes[i].FarField[ii].probe.thetastart = atan2(sqrt(pow(cos(THETASTART), 2.0) + pow(sin(PHISTART), 2) * pow(sin(THETASTART), 2)), cos(PHISTART) * sin(THETASTART));
                   this.oldSONDA.probes[i].FarField[ii].probe.phistart = atan2(cos(THETASTART), sin(PHISTART) * sin(THETASTART));
                   this.oldSONDA.probes[i].FarField[ii].probe.thetastop = atan2(sqrt(pow(cos(THETASTOP), 2.0) + pow(sin(PHISTOP), 2) * pow(sin(THETASTOP), 2)), cos(PHISTOP) * sin(THETASTOP));
                   this.oldSONDA.probes[i].FarField[ii].probe.phistop = atan2(cos(THETASTOP), sin(PHISTOP) * sin(THETASTOP));
            }
            tama3 = (this.oldSONDA.probes[i].FarField[ii].probe.n_cord);
            for (iii = 1; iii <= tama3; iii++) {
              // ROTATE MPI
              if (MPIDIR == 2) {
                  iox = this.oldSONDA.probes[i].FarField[ii].probe.i[iii];
                  ioy = this.oldSONDA.probes[i].FarField[ii].probe.j[iii];
                  ioz = this.oldSONDA.probes[i].FarField[ii].probe.K[iii];

                 this.oldSONDA.probes[i].FarField[ii].probe.i[iii] = ioz;
                 this.oldSONDA.probes[i].FarField[ii].probe.j[iii] = iox;
                 this.oldSONDA.probes[i].FarField[ii].probe.K[iii] = ioy;
              } else if (MPIDIR == 1) {
                 iox = this.oldSONDA.probes[i].FarField[ii].probe.i[iii];
                 ioy = this.oldSONDA.probes[i].FarField[ii].probe.j[iii];
                 ioz = this.oldSONDA.probes[i].FarField[ii].probe.K[iii];

                 this.oldSONDA.probes[i].FarField[ii].probe.i[iii] = ioy;
                 this.oldSONDA.probes[i].FarField[ii].probe.j[iii] = ioz;
                 this.oldSONDA.probes[i].FarField[ii].probe.K[iii] = iox;
              }
            }
            delete old_FarField;
         }
      }
      //
      for (i = 1; i <= tama; i++) {
         tama2 = (this.oldSONDA.probes[i].n_Electric);
         for (ii = 1; ii <= tama2; ii++) {
            old_Electric = new Electric_Sonda_t(this.oldSONDA.probes[i].Electric[ii]);
            tama3 = (this.oldSONDA.probes[i].Electric[ii].probe.n_cord);
            for (iii = 1; iii <= tama3; iii++) {
              // ROTATE MPI
              if (MPIDIR == 2) {
                 iox = this.oldSONDA.probes[i].Electric[ii].probe.i[iii];
                 ioy = this.oldSONDA.probes[i].Electric[ii].probe.j[iii];
                 ioz = this.oldSONDA.probes[i].Electric[ii].probe.K[iii];

                 this.oldSONDA.probes[i].Electric[ii].probe.i[iii] = ioz;
                 this.oldSONDA.probes[i].Electric[ii].probe.j[iii] = iox;
                 this.oldSONDA.probes[i].Electric[ii].probe.K[iii] = ioy;
              } else if (MPIDIR == 1) {
                 iox = this.oldSONDA.probes[i].Electric[ii].probe.i[iii];
                 ioy = this.oldSONDA.probes[i].Electric[ii].probe.j[iii];
                 ioz = this.oldSONDA.probes[i].Electric[ii].probe.K[iii];

                 this.oldSONDA.probes[i].Electric[ii].probe.i[iii] = ioy;
                 this.oldSONDA.probes[i].Electric[ii].probe.j[iii] = ioz;
                 this.oldSONDA.probes[i].Electric[ii].probe.K[iii] = iox;
              }
            }
            delete old_Electric;
         }
      }
      //
      for (i = 1; i <= tama; i++) {
         tama2 = (this.oldSONDA.probes[i].n_Magnetic);
         for (ii = 1; ii <= tama2; ii++) {
            old_Magnetic = new Magnetic_Sonda_t(this.oldSONDA.probes[i].Magnetic[ii]);
            tama3 = (this.oldSONDA.probes[i].Magnetic[ii].probe.n_cord);
            for (iii = 1; iii <= tama3; iii++) {
              // ROTATE MPI
              if (MPIDIR == 2) {
                 iox = this.oldSONDA.probes[i].Magnetic[ii].probe.i[iii];
                 ioy = this.oldSONDA.probes[i].Magnetic[ii].probe.j[iii];
                 ioz = this.oldSONDA.probes[i].Magnetic[ii].probe.K[iii];

                 this.oldSONDA.probes[i].Magnetic[ii].probe.i[iii] = ioz;
                 this.oldSONDA.probes[i].Magnetic[ii].probe.j[iii] = iox;
                 this.oldSONDA.probes[i].Magnetic[ii].probe.K[iii] = ioy;
              } else if (MPIDIR == 1) {
                 iox = this.oldSONDA.probes[i].Magnetic[ii].probe.i[iii];
                 ioy = this.oldSONDA.probes[i].Magnetic[ii].probe.j[iii];
                 ioz = this.oldSONDA.probes[i].Magnetic[ii].probe.K[iii];

                 this.oldSONDA.probes[i].Magnetic[ii].probe.i[iii] = ioy;
                 this.oldSONDA.probes[i].Magnetic[ii].probe.j[iii] = ioz;
                 this.oldSONDA.probes[i].Magnetic[ii].probe.K[iii] = iox;
              }
            }
            delete old_Magnetic;
         }
      }
      //
      return;
   }

   void rotate_generateMasSondas(Parseador_t& this, int mpidir) {
      int tama, tama2, i, ii;
      int oxi, oyi, ozi, oxe, oye, oze, oor, TXI, TYI, TZI;
      coords_t* old_MasSonda = nullptr;

      tama = this.Sonda.length;
      // tres posibilidades FarField, Electric, Magnetic
      for (i = 1; i <= tama; i++) {
         tama2 = (this.Sonda.collection[i].len_cor);
         for (ii = 1; ii <= tama2; ii++) {
              old_MasSonda = new coords_t(this.Sonda.collection[i].cordinates[ii]);
              OXI = old_MasSonda->XI;
              OXE = old_MasSonda->XE;
              OYI = old_MasSonda->YI;
              OYE = old_MasSonda->YE;
              OZI = old_MasSonda->ZI;
              OZE = old_MasSonda->ZE;
              OOR = old_MasSonda->OR;
              TXI = old_MasSonda->Xtrancos;
              TYI = old_MasSonda->Ytrancos;
              TZI = old_MasSonda->Ztrancos;
              if ((OOR != NP_COR_EX) && (OOR != NP_COR_EY) && (OOR != NP_COR_EZ) &&
              (OOR != NP_COR_HX) && (OOR != NP_COR_HY) && (OOR != NP_COR_HZ)) return;
              // LAS IW Y LAS VG NO SE ROTAN
              if (MPIDIR == 2) {
                 this.Sonda.collection[i].cordinates[ii].XI = OZI;
                 this.Sonda.collection[i].cordinates[ii].XE = OZE;
                 this.Sonda.collection[i].cordinates[ii].Xtrancos = TZI;

                 this.Sonda.collection[i].cordinates[ii].YI = OXI;
                 this.Sonda.collection[i].cordinates[ii].YE = OXE;
                 this.Sonda.collection[i].cordinates[ii].Ytrancos = TXI;

                 this.Sonda.collection[i].cordinates[ii].ZI = OYI;
                 this.Sonda.collection[i].cordinates[ii].ZE = OYE;
                 this.Sonda.collection[i].cordinates[ii].Ztrancos = TYI;

                 if (OOR == NP_COR_EX) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_EY;
                 if (OOR == NP_COR_EY) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_EZ;
                 if (OOR == NP_COR_EZ) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_EX;

                 if (OOR == NP_COR_hX) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_HY;
                 if (OOR == NP_COR_hY) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_HZ;
                 if (OOR == NP_COR_hZ) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_HX;
              } else if (MPIDIR == 1) {
                 this.Sonda.collection[i].cordinates[ii].XI = OYI;
                 this.Sonda.collection[i].cordinates[ii].XE = OYE;
                 this.Sonda.collection[i].cordinates[ii].Xtrancos = TYI;

                 this.Sonda.collection[i].cordinates[ii].YI = OZI;
                 this.Sonda.collection[i].cordinates[ii].YE = OZE;
                 this.Sonda.collection[i].cordinates[ii].Ytrancos = TZI;

                 this.Sonda.collection[i].cordinates[ii].ZI = OXI;
                 this.Sonda.collection[i].cordinates[ii].ZE = OXE;
                 this.Sonda.collection[i].cordinates[ii].Ztrancos = TXI;

                 if (OOR == NP_COR_EX) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_EZ;
                 if (OOR == NP_COR_EY) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_EX;
                 if (OOR == NP_COR_EZ) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_EY;

                 if (OOR == NP_COR_HX) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_hZ;
                 if (OOR == NP_COR_HY) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_hX;
                 if (OOR == NP_COR_HZ) this.Sonda.collection[i].cordinates[ii].OR = NP_COR_hY;
              }
              delete old_MasSonda;
         }
      }
      return;
   }

   void rotate_generateBloqueProbes(Parseador_t& this, int mpidir) {
      int tama, i;
      int oxi, oyi, ozi, oxe, oye, oze;
      BloqueProbe_t* old_BloqueProbe = nullptr;

      tama = this.BloquePRB.N_BP;
      for (i = 1; i <= tama; i++) {
          old_BloqueProbe = new BloqueProbe_t(this.BloquePRB.BP[i]);
          // MPI ROTATE Bloque CURRENT
          if (MPIDIR == 2) {
             OXI = old_BloqueProbe->i1;
             OXE = old_BloqueProbe->i2;
             OYI = old_BloqueProbe->j1;
             OYE = old_BloqueProbe->j2;
             OZI = old_BloqueProbe->k1;
             OZE = old_BloqueProbe->k2;
             //
             this.BloquePRB.BP[i].i1 = OZI;
             this.BloquePRB.BP[i].i2 = OZE;
             this.BloquePRB.BP[i].j1 = OXI;
             this.BloquePRB.BP[i].j2 = OXE;
             this.BloquePRB.BP[i].k1 = OYI;
             this.BloquePRB.BP[i].k2 = OYE;
             switch (this.BloquePRB.BP[i].nml) {
              case iEx:
                this.BloquePRB.BP[i].nml = iEy;
                break;
              case iEy:
                this.BloquePRB.BP[i].nml = iEz;
                break;
              case iEz:
                this.BloquePRB.BP[i].nml = iEx;
                break;
              default:
                break;
             }
          } else if (MPIDIR == 1) {
             OXI = old_BloqueProbe->i1;
             OXE = old_BloqueProbe->i2;
             OYI = old_BloqueProbe->j1;
             OYE = old_BloqueProbe->j2;
             OZI = old_BloqueProbe->k1;
             OZE = old_BloqueProbe->k2;
             //
             this.BloquePRB.BP[i].i1 = OYI;
             this.BloquePRB.BP[i].i2 = OYE;
             this.BloquePRB.BP[i].j1 = OZI;
             this.BloquePRB.BP[i].j2 = OZE;
             this.BloquePRB.BP[i].k1 = OXI;
             this.BloquePRB.BP[i].k2 = OXE;
             switch (this.BloquePRB.BP[i].nml) {
              case iEx:
                this.BloquePRB.BP[i].nml = iEz;
                break;
              case iEy:
                this.BloquePRB.BP[i].nml = iEx;
                break;
              case iEz:
                this.BloquePRB.BP[i].nml = iEy;
                break;
              default:
                break;
             }
          }
          delete old_BloqueProbe;
      }

      // FIN ROTATE
      return;
   }

   void rotate_generateVolumicProbes(Parseador_t& this, int mpidir) {
      int tama, tama2, i, ii;
      int oxi, oyi, ozi, oxe, oye, oze, oor, TXI, TYI, TZI;
      coords_t* old_Coordinates = nullptr;

      tama = this.VolPrb.length;
      // tres posibilidades FarField, Electric, Magnetic
      for (i = 1; i <= tama; i++) {
         tama2 = (this.VolPrb.collection[i].len_cor);
         for (ii = 1; ii <= tama2; ii++) {
              old_Coordinates = new coords_t(this.VolPrb.collection[i].cordinates[ii]);
              OXI = old_Coordinates->XI;
              OXE = old_Coordinates->XE;
              TXI = old_Coordinates->Xtrancos;
              OYI = old_Coordinates->YI;
              OYE = old_Coordinates->YE;
              TYI = old_Coordinates->Ytrancos;
              OZI = old_Coordinates->ZI;
              OZE = old_Coordinates->ZE;
              TZI = old_Coordinates->Ztrancos;
              OOR = old_Coordinates->OR;

              OXI = old_Coordinates->XI;
              OXE = old_Coordinates->XE;
              OYI = old_Coordinates->YI;
              OYE = old_Coordinates->YE;
              OZI = old_Coordinates->ZI;
              OZE = old_Coordinates->ZE;
              OOR = old_Coordinates->OR;
              // TRANCOS
              TXI = old_Coordinates->XTRANCOS;
              TYI = old_Coordinates->YTRANCOS;
              TZI = old_Coordinates->ZTRANCOS;
              //
              if ((OOR != iExC) && (OOR != iEyC) && (OOR != iEzC) &&
              (OOR != iHxC) && (OOR != iHyC) && (OOR != iHzC) &&
              (OOR != iCurX) && (OOR != iCurY) && (OOR != iCurZ) &&
               (OOR != iMEC) && (OOR != iMHC) && (OOR != iCur)) return;
              // LAS IW Y LAS VG NO SE ROTAN.
              // las imec, imhc e icur no le afecta el oor
              if (MPIDIR == 2) {
                 this.VolPrb.collection[i].cordinates[ii].XI = OZI;
                 this.VolPrb.collection[i].cordinates[ii].XE = OZE;
                 this.VolPrb.collection[i].cordinates[ii].YI = OXI;
                 this.VolPrb.collection[i].cordinates[ii].YE = OXE;
                 this.VolPrb.collection[i].cordinates[ii].ZI = OYI;
                 this.VolPrb.collection[i].cordinates[ii].ZE = OYE;
        //
                 this.VolPrb.collection[i].cordinates[ii].XTRANCOS = TZI;
                 this.VolPrb.collection[i].cordinates[ii].YTRANCOS = TXI;
                 this.VolPrb.collection[i].cordinates[ii].ZTRANCOS = TYI;
                 //
                 if (OOR == iEXc) this.VolPrb.collection[i].cordinates[ii].OR = iEYc;
                 if (OOR == iEYc) this.VolPrb.collection[i].cordinates[ii].OR = iEZc;
                 if (OOR == iEZc) this.VolPrb.collection[i].cordinates[ii].OR = iEXc;

                 if (OOR == ihXc) this.VolPrb.collection[i].cordinates[ii].OR = iHYc;
                 if (OOR == ihYc) this.VolPrb.collection[i].cordinates[ii].OR = iHZc;
                 if (OOR == ihZc) this.VolPrb.collection[i].cordinates[ii].OR = iHXc;

                 if (OOR == iCurX) this.VolPrb.collection[i].cordinates[ii].OR = iCurY;
                 if (OOR == iCurY) this.VolPrb.collection[i].cordinates[ii].OR = iCurZ;
                 if (OOR == iCurZ) this.VolPrb.collection[i].cordinates[ii].OR = iCurX;
              } else if (MPIDIR == 1) {
                 this.VolPrb.collection[i].cordinates[ii].XI = OYI;
                 this.VolPrb.collection[i].cordinates[ii].XE = OYE;
                 this.VolPrb.collection[i].cordinates[ii].YI = OZI;
                 this.VolPrb.collection[i].cordinates[ii].YE = OZE;
                 this.VolPrb.collection[i].cordinates[ii].ZI = OXI;
                 this.VolPrb.collection[i].cordinates[ii].ZE = OXE;
        //
                 this.VolPrb.collection[i].cordinates[ii].XTRANCOS = TYI;
                 this.VolPrb.collection[i].cordinates[ii].YTRANCOS = TZI;
                 this.VolPrb.collection[i].cordinates[ii].ZTRANCOS = TXI;
                 //
                 if (OOR == iEXc) this.VolPrb.collection[i].cordinates[ii].OR = iEZc;
                 if (OOR == iEYc) this.VolPrb.collection[i].cordinates[ii].OR = iEXc;
                 if (OOR == iEZc) this.VolPrb.collection[i].cordinates[ii].OR = iEYc;

                 if (OOR == iHXc) this.VolPrb.collection[i].cordinates[ii].OR = ihZc;
                 if (OOR == iHYc) this.VolPrb.collection[i].cordinates[ii].OR = ihXc;
                 if (OOR == iHZc) this.VolPrb.collection[i].cordinates[ii].OR = ihYc;

                 if (OOR == iCurX) this.VolPrb.collection[i].cordinates[ii].OR = iCurZ;
                 if (OOR == iCurY) this.VolPrb.collection[i].cordinates[ii].OR = iCurX;
                 if (OOR == iCurZ) this.VolPrb.collection[i].cordinates[ii].OR = iCurY;
              }
              delete old_Coordinates;
         }
      }

      return;
   }

   void ROTATEMPI(int mpidir, coords_t& COORDEN) {
      int OXI, OXE, OYI, OYE, OZI, OZE, OOR, TXI, TYI, TZI;
      OXI = COORDEN.XI;
      OXE = COORDEN.XE;
      TXI = COORDEN.Xtrancos;
      OYI = COORDEN.YI;
      OYE = COORDEN.YE;
      TYI = COORDEN.Ytrancos;
      OZI = COORDEN.ZI;
      OZE = COORDEN.ZE;
      TZI = COORDEN.Ztrancos;
      OOR = COORDEN.OR;
      if (MPIDIR == 2) {
         COORDEN.XI = OZI;
         COORDEN.XE = OZE;
         COORDEN.Xtrancos = TZI;
         COORDEN.YI = OXI;
         COORDEN.YE = OXE;
         COORDEN.Ytrancos = TXI;
         COORDEN.ZI = OYI;
         COORDEN.ZE = OYE;
         COORDEN.Ztrancos = TYI;
         if (OOR == iEx) COORDEN.OR = iEy;
         if (OOR == -iEx) COORDEN.OR = -iEy;
         if (OOR == iEy) COORDEN.OR = iEz;
         if (OOR == -iEy) COORDEN.OR = -iEz;
         if (OOR == iEz) COORDEN.OR = iEx;
         if (OOR == -iEz) COORDEN.OR = -iEx;
      } else if (MPIDIR == 1) {
         COORDEN.XI = OYI;
         COORDEN.XE = OYE;
         COORDEN.Xtrancos = TYI;
         COORDEN.YI = OZI;
         COORDEN.YE = OZE;
         COORDEN.Ytrancos = TZI;
         COORDEN.ZI = OXI;
         COORDEN.ZE = OXE;
         COORDEN.Ztrancos = TXI;
         if (OOR == iEx) COORDEN.OR = iEz;
         if (OOR == -iEx) COORDEN.OR = -iEz;
         if (OOR == iEy) COORDEN.OR = iEx;
         if (OOR == -iEy) COORDEN.OR = -iEx;
         if (OOR == iEz) COORDEN.OR = iEy;
         if (OOR == -iEz) COORDEN.OR = -iEy;
      }
      return;
   }

   void ROTATEMPI_SCALED(int mpidir, coords_scaled_t& COORDEN) {
      int OXI, OXE, OYI, OYE, OZI, OZE, oor;
      double OXC, OYC, OZC;
      OXI = COORDEN.XI;
      OXE = COORDEN.XE;
      OYI = COORDEN.YI;
      OYE = COORDEN.YE;
      OZI = COORDEN.ZI;
      OZE = COORDEN.ZE;
      OXC = COORDEN.XC;
      OYC = COORDEN.YC;
      OZC = COORDEN.ZC;
      OOr = COORDEN.or;
      if (MPIDIR == 2) {
         COORDEN.XI = OZI;
         COORDEN.XE = OZE;
         COORDEN.YI = OXI;
         COORDEN.YE = OXE;
         COORDEN.ZI = OYI;
         COORDEN.ZE = OYE;
         COORDEN.XC = OzC;
         COORDEN.YC = OxC;
         COORDEN.ZC = OyC;
         if (OOR == iEx) COORDEN.OR = iEy;
         if (OOR == -iEx) COORDEN.OR = -iEy;
         if (OOR == iEy) COORDEN.OR = iEz;
         if (OOR == -iEy) COORDEN.OR = -iEz;

if (OOR == iEz) COORDEN.OR = iEx;
        if (OOR == -iEz) COORDEN.OR = -iEx;

    } else if (MPIDIR == 1) {
        COORDEN.XI = OYI;
        COORDEN.XE = OYE;
        COORDEN.YI = OZI;
        COORDEN.YE = OZE;
        COORDEN.ZI = OXI;
        COORDEN.ZE = OXE;
        COORDEN.XC = OyC;
        COORDEN.YC = OzC;
        COORDEN.ZC = OxC;
        if (OOR == iEx) COORDEN.OR = iEz;
        if (OOR == -iEx) COORDEN.OR = -iEz;
        if (OOR == iEy) COORDEN.OR = iEx;
        if (OOR == -iEy) COORDEN.OR = -iEx;
        if (OOR == iEz) COORDEN.OR = iEy;
        if (OOR == -iEz) COORDEN.OR = -iEy;
    }
    return;
}

void rotate_freq_depend_material_properties(int mpidir, FreqDepenMaterial_t& freqDepMat) {
    std::complex<double>* po11 = nullptr;
    std::complex<double>* po12 = nullptr;
    std::complex<double>* po13 = nullptr;
    std::complex<double>* po22 = nullptr;
    std::complex<double>* po23 = nullptr;
    std::complex<double>* po33 = nullptr;
    double ro11, ro12, ro13, ro22, ro23, ro33;
    int io11, io12, io13, io22, io23, io33;

    if (mpidir == 2) {
        // Rotate a matrix
        po11 = &freqDepMat.a11;
        po12 = &freqDepMat.a12;
        po13 = &freqDepMat.a13;
        po22 = &freqDepMat.a22;
        po23 = &freqDepMat.a23;
        po33 = &freqDepMat.a33;

        freqDepMat.a11 = *po33;
        freqDepMat.a12 = *po23;
        freqDepMat.a13 = *po12;
        freqDepMat.a22 = *po11;
        freqDepMat.a23 = *po13;
        freqDepMat.a33 = *po22;

        // Rotate am matrix
        po11 = &freqDepMat.am11;
        po12 = &freqDepMat.am12;
        po13 = &freqDepMat.am13;
        po22 = &freqDepMat.am22;
        po23 = &freqDepMat.am23;
        po33 = &freqDepMat.am33;

        freqDepMat.am11 = *po33;
        freqDepMat.am12 = *po23;
        freqDepMat.am13 = *po12;
        freqDepMat.am22 = *po11;
        freqDepMat.am23 = *po13;
        freqDepMat.am33 = *po22;

        // Rotate b matrix
        po11 = &freqDepMat.b11;
        po12 = &freqDepMat.b12;
        po13 = &freqDepMat.b13;
        po22 = &freqDepMat.b22;
        po23 = &freqDepMat.b23;
        po33 = &freqDepMat.b33;

        freqDepMat.b11 = *po33;
        freqDepMat.b12 = *po23;
        freqDepMat.b13 = *po12;
        freqDepMat.b22 = *po11;
        freqDepMat.b23 = *po13;
        freqDepMat.b33 = *po22;

        // Rotate bm matrix
        po11 = &freqDepMat.bm11;
        po12 = &freqDepMat.bm12;
        po13 = &freqDepMat.bm13;
        po22 = &freqDepMat.bm22;
        po23 = &freqDepMat.bm23;
        po33 = &freqDepMat.bm33;

        freqDepMat.bm11 = *po33;
        freqDepMat.bm12 = *po23;
        freqDepMat.bm13 = *po12;
        freqDepMat.bm22 = *po11;
        freqDepMat.bm23 = *po13;
        freqDepMat.bm33 = *po22;

        // Rotate eps matrix
        ro11 = freqDepMat.eps11;
        ro12 = freqDepMat.eps12;
        ro13 = freqDepMat.eps13;
        ro22 = freqDepMat.eps22;
        ro23 = freqDepMat.eps23;
        ro33 = freqDepMat.eps33;

        freqDepMat.eps11 = ro33;
        freqDepMat.eps12 = ro23;
        freqDepMat.eps13 = ro12;
        freqDepMat.eps22 = ro11;
        freqDepMat.eps23 = ro13;
        freqDepMat.eps33 = ro22;

        // Rotate mu matrix
        ro11 = freqDepMat.mu11;
        ro12 = freqDepMat.mu12;
        ro13 = freqDepMat.mu13;
        ro22 = freqDepMat.mu22;
        ro23 = freqDepMat.mu23;
        ro33 = freqDepMat.mu33;

        freqDepMat.mu11 = ro33;
        freqDepMat.mu12 = ro23;
        freqDepMat.mu13 = ro12;
        freqDepMat.mu22 = ro11;
        freqDepMat.mu23 = ro13;
        freqDepMat.mu33 = ro22;

        // Rotate sigma matrix
        ro11 = freqDepMat.sigma11;
        ro12 = freqDepMat.sigma12;
        ro13 = freqDepMat.sigma13;
        ro22 = freqDepMat.sigma22;
        ro23 = freqDepMat.sigma23;
        ro33 = freqDepMat.sigma33;

        freqDepMat.sigma11 = ro33;
        freqDepMat.sigma12 = ro23;
        freqDepMat.sigma13 = ro12;
        freqDepMat.sigma22 = ro11;
        freqDepMat.sigma23 = ro13;
        freqDepMat.sigma33 = ro22;

        // Rotate sigmam matrix
        ro11 = freqDepMat.sigmam11;
        ro12 = freqDepMat.sigmam12;
        ro13 = freqDepMat.sigmam13;
        ro22 = freqDepMat.sigmam22;
        ro23 = freqDepMat.sigmam23;
        ro33 = freqDepMat.sigmam33;

        freqDepMat.sigmam11 = ro33;
        freqDepMat.sigmam12 = ro23;
        freqDepMat.sigmam13 = ro12;
        freqDepMat.sigmam22 = ro11;
        freqDepMat.sigmam23 = ro13;
        freqDepMat.sigmam33 = ro22;

        // Rotate K matrix
        io11 = freqDepMat.k11;
        io12 = freqDepMat.k12;
        io13 = freqDepMat.k13;
        io22 = freqDepMat.k22;
        io23 = freqDepMat.k23;
        io33 = freqDepMat.k33;

        freqDepMat.k11 = io33;
        freqDepMat.k12 = io23;
        freqDepMat.k13 = io12;
        freqDepMat.k22 = io11;
        freqDepMat.k23 = io13;
        freqDepMat.k33 = io22;

        // Rotate Km matrix
        io11 = freqDepMat.km11;
        io12 = freqDepMat.km12;
        io13 = freqDepMat.km13;
        io22 = freqDepMat.km22;
        io23 = freqDepMat.km23;
        io33 = freqDepMat.km33;

        freqDepMat.km11 = io33;
        freqDepMat.km12 = io23;
        freqDepMat.km13 = io12;
        freqDepMat.km22 = io11;
        freqDepMat.km23 = io13;
        freqDepMat.km33 = io22;
    }

    if (mpidir == 1) {
        // Rotate a matrix
        po11 = &freqDepMat.a11;
        po12 = &freqDepMat.a12;
        po13 = &freqDepMat.a13;
        po22 = &freqDepMat.a22;
        po23 = &freqDepMat.a23;
        po33 = &freqDepMat.a33;

        freqDepMat.a11 = *po22;
        freqDepMat.a12 = *po13;
        freqDepMat.a13 = *po23;
        freqDepMat.a22 = *po33;
        freqDepMat.a23 = *po12;
        freqDepMat.a33 = *po11;

        // Rotate am matrix
        po11 = &freqDepMat.am11;
        po12 = &freqDepMat.am12;
        po13 = &freqDepMat.am13;
        po22 = &freqDepMat.am22;
        po23 = &freqDepMat.am23;
        po33 = &freqDepMat.am33;

        freqDepMat.am11 = *po22;
        freqDepMat.am12 = *po13;
        freqDepMat.am13 = *po23;
        freqDepMat.am22 = *po33;
        freqDepMat.am23 = *po12;
        freqDepMat.am33 = *po11;

        // Rotate b matrix
        po11 = &freqDepMat.b11;
        po12 = &freqDepMat.b12;
        po13 = &freqDepMat.b13;
        po22 = &freqDepMat.b22;
        po23 = &freqDepMat.b23;
        po33 = &freqDepMat.b33;

        freqDepMat.b11 = *po22;
        freqDepMat.b12 = *po13;
        freqDepMat.b13 = *po23;
        freqDepMat.b22 = *po33;
        freqDepMat.b23 = *po12;
        freqDepMat.b33 = *po11;

        // Rotate bm matrix
        po11 = &freqDepMat.bm11;
        po12 = &freqDepMat.bm12;
        po13 = &freqDepMat.bm13;
        po22 = &freqDepMat.bm22;
        po23 = &freqDepMat.bm23;
        po33 = &freqDepMat.bm33;

        freqDepMat.bm11 = *po22;
        freqDepMat.bm12 = *po13;
        freqDepMat.bm13 = *po23;
        freqDepMat.bm22 = *po33;
        freqDepMat.bm23 = *po12;
        freqDepMat.bm33 = *po11;

        // Rotate eps matrix
        ro11 = freqDepMat.eps11;
        ro12 = freqDepMat.eps12;
        ro13 = freqDepMat.eps13;
        ro22 = freqDepMat.eps22;
        ro23 = freqDepMat.eps23;
        ro33 = freqDepMat.eps33;

        freqDepMat.eps11 = ro22;
        freqDepMat.eps12 = ro13;
        freqDepMat.eps13 = ro23;
        freqDepMat.eps22 = ro33;
        freqDepMat.eps23 = ro12;
        freqDepMat.eps33 = ro11;

        // Rotate mu matrix
        ro11 = freqDepMat.mu11;
        ro12 = freqDepMat.mu12;
        ro13 = freqDepMat.mu13;
        ro22 = freqDepMat.mu22;
        ro23 = freqDepMat.mu23;
        ro33 = freqDepMat.mu33;

        freqDepMat.mu11 = ro22;
        freqDepMat.mu12 = ro13;
        freqDepMat.mu13 = ro23;
        freqDepMat.mu22 = ro33;
        freqDepMat.mu23 = ro12;
        freqDepMat.mu33 = ro11;

        // Rotate sigma matrix
        ro11 = freqDepMat.sigma11;
        ro12 = freqDepMat.sigma12;
        ro13 = freqDepMat.sigma13;
        ro22 = freqDepMat.sigma22;
        ro23 = freqDepMat.sigma23;
        ro33 = freqDepMat.sigma33;

        freqDepMat.sigma11 = ro22;
        freqDepMat.sigma12 = ro13;
        freqDepMat.sigma13 = ro23;
        freqDepMat.sigma22 = ro33;
        freqDepMat.sigma23 = ro12;
        freqDepMat.sigma33 = ro11;

        // Rotate sigmam matrix
        ro11 = freqDepMat.sigmam11;
        ro12 = freqDepMat.sigmam12;
        ro13 = freqDepMat.sigmam13;
        ro22 = freqDepMat.sigmam22;
        ro23 = freqDepMat.sigmam23;
        ro33 = freqDepMat.sigmam33;

        freqDepMat.sigmam11 = ro22;
        freqDepMat.sigmam12 = ro13;
        freqDepMat.sigmam13 = ro23;
        freqDepMat.sigmam22 = ro33;
        freqDepMat.sigmam23 = ro12;
        freqDepMat.sigmam33 = ro11;

        // Rotate K matrix
        io11 = freqDepMat.k11;
        io12 = freqDepMat.k12;
        io13 = freqDepMat.k13;
        io22 = freqDepMat.k22;
        io23 = freqDepMat.k23;
        io33 = freqDepMat.k33;

        freqDepMat.k11 = io22;
        freqDepMat.k12 = io13;
        freqDepMat.k13 = io23;
        freqDepMat.k22 = io33;
        freqDepMat.k23 = io12;
        freqDepMat.k33 = io11;

        // Rotate Km matrix
        io11 = freqDepMat.km11;
        io12 = freqDepMat.km12;
        io13 = freqDepMat.km13;
        io22 = freqDepMat.km22;
        io23 = freqDepMat.km23;
        io33 = freqDepMat.km33;

        freqDepMat.km11 = io22;
        freqDepMat.km12 = io13;
        freqDepMat.km13 = io23;
        freqDepMat.km22 = io33;
        freqDepMat.km23 = io12;
        freqDepMat.km33 = io11;
    }
}