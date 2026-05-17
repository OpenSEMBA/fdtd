```cpp
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <memory>
#include <cassert>

// Forward declarations and includes for external types used in the module
// These would typically come from FDETYPES_m, Report_m, etc.
// Assuming standard definitions based on context:

// Placeholder for external types/constants that are not defined in this snippet
// In a real translation, these would be included from their respective headers
namespace FDETYPES_m {
    using RKIND = double;
    using INTEGERSIZEOFMEDIAMATRICES = int;
    
    // Placeholder for SGGFDTDINFO_t structure
    struct SGGFDTDINFO_t {
        struct {
            struct { int XI, XE, YI, YE, ZI, ZE; } iHx, iHy, iHz, iEx, iEy, iEz;
        } ALLOC;
        
        struct {
            struct { int XI, XE, YI, YE, ZI, ZE; } iHx, iHy, iHz, iEx, iEy, iEz;
        } Sweep;
        
        struct {
            int NumLayers[4][3]; // Assuming dimensions based on usage o=1..3, p=1..2
            int orden[4][3];
            double CoeffReflPML[4][3];
        } PML;
        
        struct {
            bool IsBackPML, IsFrontPML, IsLeftPML, IsRightPML, IsUpPML, IsDownPML;
        } Border;
        
        double dt;
    };
    
    // Placeholder for sim_control_t
    struct sim_control_t {
        double alphamaxpar;
        int alphaOrden;
        double kappamaxpar;
        bool resume;
    };
    
    // Placeholder for bounds_t
    struct bounds_t {
        struct { int NX, NY, NZ, XI, XE, YI, YE, ZI, ZE; } sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz, Ex, Ey, Ez, Hx, Hy, Hz;
    };
    
    // Placeholder for limit_t
    struct limit_t {
        int XI, XE, YI, YE, ZI, ZE;
    };
    
    // Constants
    constexpr int BUFSIZE = 256;
    constexpr int iEx = 1, iEy = 2, iEz = 3, iHx = 4, iHy = 5, iHz = 6;
    constexpr int Down = 1, Up = 2, Left = 3, Right = 4, Back = 5, Front = 6;
    
    // External functions/variables placeholders
    void print11(int layoutnumber, const std::string& message);
    extern std::string SEPARADOR;
}

using namespace FDETYPES_m;

namespace BORDERS_CPML_m {

    // Constants
    constexpr double StaticFrequency = 1.0e14;

    // Derived types converted to structs
    struct xyzlimit_var_t {
        int XI[6], XE[6], YI[6], YE[6], ZI[6], ZE[6]; // 1-based indexing mapped to 0-based array with offset or direct mapping
        // Fortran: dimension(1:6). We'll use 0-based vector or array and adjust access.
        // To preserve logic easily, let's use std::vector<int> of size 7, ignoring index 0.
        std::vector<int> XI_vec, XE_vec, YI_vec, YE_vec, ZI_vec, ZE_vec;
        
        xyzlimit_var_t() : XI_vec(7), XE_vec(7), YI_vec(7), YE_vec(7), ZI_vec(7), ZE_vec(7) {}
        
        int& XI(int idx) { return XI_vec[idx]; }
        int& XE(int idx) { return XE_vec[idx]; }
        int& YI(int idx) { return YI_vec[idx]; }
        int& YE(int idx) { return YE_vec[idx]; }
        int& ZI(int idx) { return ZI_vec[idx]; }
        int& ZE(int idx) { return ZE_vec[idx]; }
        
        const int& XI(int idx) const { return XI_vec[idx]; }
        const int& XE(int idx) const { return XE_vec[idx]; }
        const int& YI(int idx) const { return YI_vec[idx]; }
        const int& YE(int idx) const { return YE_vec[idx]; }
        const int& ZI(int idx) const { return ZI_vec[idx]; }
        const int& ZE(int idx) const { return ZE_vec[idx]; }
    };

    struct LR_t {
        // 3D arrays. Using flattened 1D vector for performance and ease of management, 
        // or 3D vector. Given the indexing, a 3D vector or a helper class is best.
        // For simplicity and direct translation of indexing, we'll use a struct with vectors.
        // However, to match Fortran's dynamic allocation, std::vector is appropriate.
        // We will use a helper to manage 3D indexing.
        
        struct Array3D {
            std::vector<double> data;
            int nx, ny, nz;
            int ox, oy, oz; // offsets for 1-based or arbitrary start indices
            
            Array3D() : nx(0), ny(0), nz(0), ox(0), oy(0), oz(0) {}
            
            void allocate(int x1, int x2, int y1, int y2, int z1, int z2) {
                nx = x2 - x1 + 1;
                ny = y2 - y1 + 1;
                nz = z2 - z1 + 1;
                ox = x1; oy = y1; oz = z1;
                data.resize(nx * ny * nz, 0.0);
            }
            
            double& operator()(int i, int j, int k) {
                return data[(k - oz) * nx * ny + (j - oy) * nx + (i - ox)];
            }
            
            const double& operator()(int i, int j, int k) const {
                return data[(k - oz) * nx * ny + (j - oy) * nx + (i - ox)];
            }
            
            void fill(double val) {
                std::fill(data.begin(), data.end(), val);
            }
        };

        Array3D Psi_Exy, Psi_Ezy, Psi_Hxy, Psi_Hzy;
        Array3D Psi_Exyvac, Psi_Ezyvac, Psi_Hxyvac, Psi_Hzyvac;
    };

    struct DU_t {
        struct Array3D {
            std::vector<double> data;
            int nx, ny, nz;
            int ox, oy, oz;
            
            Array3D() : nx(0), ny(0), nz(0), ox(0), oy(0), oz(0) {}
            
            void allocate(int x1, int x2, int y1, int y2, int z1, int z2) {
                nx = x2 - x1 + 1;
                ny = y2 - y1 + 1;
                nz = z2 - z1 + 1;
                ox = x1; oy = y1; oz = z1;
                data.resize(nx * ny * nz, 0.0);
            }
            
            double& operator()(int i, int j, int k) {
                return data[(k - oz) * nx * ny + (j - oy) * nx + (i - ox)];
            }
            
            const double& operator()(int i, int j, int k) const {
                return data[(k - oz) * nx * ny + (j - oy) * nx + (i - ox)];
            }
            
            void fill(double val) {
                std::fill(data.begin(), data.end(), val);
            }
        };

        Array3D Psi_Eyz, Psi_Exz, Psi_Hyz, Psi_Hxz;
        Array3D Psi_Eyzvac, Psi_Exzvac, Psi_Hyzvac, Psi_Hxzvac;
    };

    struct BF_t {
        struct Array3D {
            std::vector<double> data;
            int nx, ny, nz;
            int ox, oy, oz;
            
            Array3D() : nx(0), ny(0), nz(0), ox(0), oy(0), oz(0) {}
            
            void allocate(int x1, int x2, int y1, int y2, int z1, int z2) {
                nx = x2 - x1 + 1;
                ny = y2 - y1 + 1;
                nz = z2 - z1 + 1;
                ox = x1; oy = y1; oz = z1;
                data.resize(nx * ny * nz, 0.0);
            }
            
            double& operator()(int i, int j, int k) {
                return data[(k - oz) * nx * ny + (j - oy) * nx + (i - ox)];
            }
            
            const double& operator()(int i, int j, int k) const {
                return data[(k - oz) * nx * ny + (j - oy) * nx + (i - ox)];
            }
            
            void fill(double val) {
                std::fill(data.begin(), data.end(), val);
            }
        };

        Array3D Psi_Ezx, Psi_Eyx, Psi_Hzx, Psi_Hyx;
        Array3D Psi_Ezxvac, Psi_Eyxvac, Psi_Hzxvac, Psi_Hyxvac;
    };

    // Global variables from the module
    std::vector<LR_t> regLR(7); // left:right, assuming 1-based 6 regions + padding or just 6. Fortran says left:right. Usually 1..6.
    std::vector<DU_t> regDU(7); // down:up
    std::vector<BF_t> regBF(7); // back:front
    
    std::vector<double> sig_max; // 1:3, 1:2
    std::vector<double> aPar_max; // 1:3, 1:2
    std::vector<double> kPar_max; // 1:3, 1:2
    
    std::vector<double> P_ce_x, P_ce_y, P_ce_z;
    std::vector<double> P_be_x, P_be_y, P_be_z;
    std::vector<double> P_cm_x, P_cm_y, P_cm_z;
    std::vector<double> P_bm_x, P_bm_y, P_bm_z;
    
    std::vector<double> ce_x, ce_y, ce_z;
    std::vector<double> cm_x, cm_y, cm_z;
    std::vector<double> Ice_x, Ice_y, Ice_z;
    std::vector<double> Icm_x, Icm_y, Icm_z;
    
    double zvac = 0.0;
    double eps0 = 0.0;
    double mu0 = 0.0;
    
    double alphamaxpar = 0.0;
    int alphaOrden = 0;
    double kappamaxpar = 0.0;
    
    std::vector<limit_t> SINPML_fullsize(7); // 1:6
    
    std::vector<double> dxe, dye, dze, dxh, dyh, dzh;

    // Helper to access sig_max, aPar_max, kPar_max with 1-based indexing
    inline double& sig_max(int i, int j) { return sig_max[(i-1)*2 + (j-1)]; }
    inline double& aPar_max(int i, int j) { return aPar_max[(i-1)*2 + (j-1)]; }
    inline double& kPar_max(int i, int j) { return kPar_max[(i-1)*2 + (j-1)]; }
    inline const double& sig_max(int i, int j) const { return sig_max[(i-1)*2 + (j-1)]; }
    inline const double& aPar_max(int i, int j) const { return aPar_max[(i-1)*2 + (j-1)]; }
    inline const double& kPar_max(int i, int j) const { return kPar_max[(i-1)*2 + (j-1)]; }

    void InitCPMLBorders(const SGGFDTDINFO_t& sgg, const limit_t temp_SINPML_Fullsize[7], bool& ThereArePMLBorders, const sim_control_t& control,
                         const std::vector<double>& temp_dxe, const std::vector<double>& temp_dye, const std::vector<double>& temp_dze,
                         const std::vector<double>& temp_dxh, const std::vector<double>& temp_dyh, const std::vector<double>& temp_dzh,
                         std::vector<double>& Idxe, std::vector<double>& Idye, std::vector<double>& Idze,
                         std::vector<double>& Idxh, std::vector<double>& Idyh, std::vector<double>& Idzh,
                         double eps00, double mu00) {
        eps0 = eps00;
        mu0 = mu00;
        zvac = std::sqrt(mu0 / eps0);
        
        for (int i = 1; i <= 6; ++i) {
            SINPML_fullsize[i] = temp_SINPML_Fullsize[i];
        }
        
        alphamaxpar = control.alphamaxpar;
        alphaOrden = control.alphaOrden;
        kappamaxpar = control.kappamaxpar;
        
        // Allocate dxe, etc.
        // Fortran arrays are passed with specific bounds. We assume the vectors passed are sized correctly and indexed 1-based relative to the bounds.
        // For simplicity in C++, we'll assume the input vectors are 0-indexed but represent the data.
        // However, to preserve logic, we need to map the Fortran bounds to C++ indices.
        // The input vectors temp_dxe etc. are passed with bounds sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE.
        // We will store them in vectors that are sized to cover the range, using an offset.
        
        auto allocate_vector_with_bounds = [](const std::vector<double>& src, int xi, int xe) -> std::vector<double> {
            std::vector<double> dest(xe - xi + 1);
            for (int i = 0; i < src.size(); ++i) {
                dest[i] = src[i];
            }
            return dest;
        };

        // Note: The actual implementation of passing arrays with bounds in C++ is complex.
        // For this translation, we assume the vectors passed are already aligned or we copy them into internal vectors with offset management.
        // To keep it simple and functional, we'll copy the data into our global vectors, assuming the caller handles the indexing correctly or we adjust.
        
        // Re-allocating global vectors based on sgg bounds
        int hxi = sgg.ALLOC.iHx.XI, hxe = sgg.ALLOC.iHx.XE;
        int hyi = sgg.ALLOC.iHy.YI, hye = sgg.ALLOC.iHy.YE;
        int hzi = sgg.ALLOC.iHz.ZI, hze = sgg.ALLOC.iHz.ZE;
        int exi = sgg.ALLOC.iEx.XI, exe = sgg.ALLOC.iEx.XE;
        int eyi = sgg.ALLOC.iEy.YI, eye = sgg.ALLOC.iEy.YE;
        int ezi = sgg.ALLOC.iEz.ZI, eze = sgg.ALLOC.iEz.ZE;

        dxe.resize(hxe - hxi + 1);
        dye.resize(hye - hyi + 1);
        dze.resize(hze - hzi + 1);
        dxh.resize(exe - exi + 1);
        dyh.resize(eye - eyi + 1);
        dzh.resize(eze - ezi + 1);

        // Copy data. Assuming temp vectors are 0-indexed and correspond to the range.
        for(int i=0; i<dxe.size(); ++i) dxe[i] = temp_dxe[i];
        for(int i=0; i<dye.size(); ++i) dye[i] = temp_dye[i];
        for(int i=0; i<dze.size(); ++i) dze[i] = temp_dze[i];
        for(int i=0; i<dxh.size(); ++i) dxh[i] = temp_dxh[i];
        for(int i=0; i<dyh.size(); ++i) dyh[i] = temp_dyh[i];
        for(int i=0; i<dzh.size(); ++i) dzh[i] = temp_dzh[i];

        ThereArePMLBorders = sgg.Border.IsBackPML || sgg.Border.IsFrontPML || sgg.Border.IsLeftPML || 
                             sgg.Border.IsRightPML || sgg.Border.IsUpPML || sgg.Border.IsDownPML;
        
        if (!ThereArePMLBorders) return;

        // Find limits
        for (int field = iEx; field <= iHz; ++field) {
            // Down
            regLR[field].Psi_Exy.allocate(0,0,0,0,0,0); // Placeholder, actual allocation happens later
            // Simplified: We need to map field to indices. 
            // In Fortran: PMLc(field)%XI(Down) ...
            // We'll use a helper to set limits.
            
            // This part is complex due to the 2D array PMLc(field, region).
            // We'll assume PMLc is a vector of xyzlimit_var_t of size 7 (1-6).
            // And the limits are stored in the xyzlimit_var_t struct.
            
            // To avoid massive code duplication in this thought block, we'll implement the logic directly.
            // Note: The Fortran code uses PMLc(field)%XI(Down). 
            // We need to store these limits. Let's add a 2D array of limits to the module globals or handle it within the structs.
            // For this translation, we'll assume the limits are stored in a global 2D array `PMLc_limits[7][7]` (field, region).
            // But the original code uses `type(xyzlimit_var_t), dimension(1:6) :: PMLc`.
            // So PMLc is an array of structs. Each struct has XI, XE etc. which are arrays of size 6 (regions).
            // So PMLc(field).XI(region) gives the limit.
            
            // We need to update the global PMLc vector.
            // Since PMLc is global, we update it here.
            
            // Down
            PMLc[field].XI(1) = sgg.Sweep.get(field).XI; // Assuming get(field) returns the sweep info for the field
            PMLc[field].XE(1) = sgg.Sweep.get(field).XE;
            PMLc[field].YI(1) = sgg.Sweep.get(field).YI;
            PMLc[field].YE(1) = sgg.Sweep.get(field).YE;
            PMLc[field].ZI(1) = sgg.Sweep.get(field).ZI;
            PMLc[field].ZE(1) = std::min(SINPML_fullsize[field].ZI - 1, sgg.Sweep.get(field).ZE);
            
            // Up
            PMLc[field].XI(2) = sgg.Sweep.get(field).XI;
            PMLc[field].XE(2) = sgg.Sweep.get(field).XE;
            PMLc[field].YI(2) = sgg.Sweep.get(field).YI;
            PMLc[field].YE(2) = sgg.Sweep.get(field).YE;
            PMLc[field].ZI(2) = std::max(SINPML_fullsize[field].ZE + 1, sgg.Sweep.get(field).ZI);
            PMLc[field].ZE(2) = sgg.Sweep.get(field).ZE;
            
            // Left
            PMLc[field].XI(3) = sgg.Sweep.get(field).XI;
            PMLc[field].XE(3) = sgg.Sweep.get(field).XE;
            PMLc[field].YI(3) = sgg.Sweep.get(field).YI;
            PMLc[field].YE(3) = std::min(SINPML_fullsize[field].YI - 1, sgg.Sweep.get(field).YE);
            PMLc[field].ZI(3) = sgg.Sweep.get(field).ZI;
            PMLc[field].ZE(3) = sgg.Sweep.get(field).ZE;
            
            // Right
            PMLc[field].XI(4) = sgg.Sweep.get(field).XI;
            PMLc[field].XE(4) = sgg.Sweep.get(field).XE;
            PMLc[field].YI(4) = std::max(SINPML_fullsize[field].YE + 1, sgg.Sweep.get(field).YI);
            PMLc[field].YE(4) = sgg.Sweep.get(field).YE;
            PMLc[field].ZI(4) = sgg.Sweep.get(field).ZI;
            PMLc[field].ZE(4) = sgg.Sweep.get(field).ZE;
            
            // Back
            PMLc[field].XI(5) = sgg.Sweep.get(field).XI;
            PMLc[field].XE(5) = std::min(SINPML_fullsize[field].XI - 1, sgg.Sweep.get(field).XE);
            PMLc[field].YI(5) = sgg.Sweep.get(field).YI;
            PMLc[field].YE(5) = sgg.Sweep.get(field).YE;
            PMLc[field].ZI(5) = sgg.Sweep.get(field).ZI;
            PMLc[field].ZE(5) = sgg.Sweep.get(field).ZE;
            
            // Front
            PMLc[field].XI(6) = std::max(SINPML_fullsize[field].XE + 1, sgg.Sweep.get(field).XI);
            PMLc[field].XE(6) = sgg.Sweep.get(field).XE;
            PMLc[field].YI(6) = sgg.Sweep.get(field).YI;
            PMLc[field].YE(6) = sgg.Sweep.get(field).YE;
            PMLc[field].ZI(6) = sgg.Sweep.get(field).ZI;
            PMLc[field].ZE(6) = sgg.Sweep.get(field).ZE;
        }

        // Allocate sig_max, aPar_max, kPar_max
        sig_max.resize(6);
        aPar_max.resize(6);
        kPar_max.resize(6);

        // Allocate P_ce_x, etc.
        P_ce_x.resize(hxe - hxi + 1);
        P_ce_y.resize(hye - hyi + 1);
        P_ce_z.resize(hze - hzi + 1);
        P_be_x.resize(hxe - hxi + 1);
        P_be_y.resize(hye - hyi + 1);
        P_be_z.resize(hze - hzi + 1);
        P_cm_x.resize(hxe - hxi + 1);
        P_cm_y.resize(hye - hyi + 1);
        P_cm_z.resize(hze - hzi + 1);
        P_bm_x.resize(hxe - hxi + 1);
        P_bm_y.resize(hye - hyi + 1);
        P_bm_z.resize(hze - hzi + 1);

        ce_x.resize(hxe - hxi + 1);
        ce_y.resize(hye - hyi + 1);
        ce_z.resize(hze - hzi + 1);
        cm_x.resize(hxe - hxi + 1);
        cm_y.resize(hye - hyi + 1);
        cm_z.resize(hze - hzi + 1);

        Ice_x.resize(hxe - hxi + 1);
        Ice_y.resize(hye - hyi + 1);
        Ice_z.resize(hze - hzi + 1);
        Icm_x.resize(hxe - hxi + 1);
        Icm_y.resize(hye - hyi + 1);
        Icm_z.resize(hze - hzi + 1);

        std::fill(P_ce_x.begin(), P_ce_x.end(), 0.0);
        std::fill(P_ce_y.begin(), P_ce_y.end(), 0.0);
        std::fill(P_ce_z.begin(), P_ce_z.end(), 0.0);
        std::fill(P_cm_x.begin(), P_cm_x.end(), 0.0);
        std::fill(P_cm_y.begin(), P_cm_y.end(), 0.0);
        std::fill(P_cm_z.begin(), P_cm_z.end(), 0.0);
        std::fill(Ice_x.begin(), Ice_x.end(), 0.0);
        std::fill(Ice_y.begin(), Ice_y.end(), 0.0);
        std::fill(Ice_z.begin(), Ice_z.end(), 0.0);
        std::fill(Icm_x.begin(), Icm_x.end(), 0.0);
        std::fill(Icm_y.begin(), Icm_y.end(), 0.0);
        std::fill(Icm_z.begin(), Icm_z.end(), 0.0);

        // Depth information matrices
        for (int i = hxi; i <= hxe; ++i) {
            int idx = i - hxi;
            if (i <= SINPML_fullsize[iHx].XI && sgg.PML.NumLayers[1][1] != 0) {
                ce_x[idx] = 1.0 * (SINPML_fullsize[iHx].XI - i) / sgg.PML.NumLayers[1][1];
                Ice_x[idx] = 1.0 * (sgg.PML.NumLayers[1][1] - (SINPML_fullsize[iHx].XI - i)) / sgg.PML.NumLayers[1][1];
            } else if (i >= SINPML_fullsize[iHx].XE && sgg.PML.NumLayers[1][2] != 0) {
                ce_x[idx] = 1.0 * (i - SINPML_fullsize[iHx].XE) / sgg.PML.NumLayers[1][2];
                Ice_x[idx] = 1.0 * (sgg.PML.NumLayers[1][2] - (i - SINPML_fullsize[iHx].XE)) / sgg.PML.NumLayers[1][2];
            } else {
                ce_x[idx] = 0.0;
                Ice_x[idx] = 0.0;
            }
        }
        for (int i = hxi; i <= hxe; ++i) {
            int idx = i - hxi;
            if (i <= SINPML_fullsize[iHx].XI - 1 && sgg.PML.NumLayers[1][1] != 0) {
                cm_x[idx] = 1.0 * (SINPML_fullsize[iHx].XI - (i + 0.5)) / sgg.PML.NumLayers[1][1];
                Icm_x[idx] = 1.0 * (sgg.PML.NumLayers[1][1] - (SINPML_fullsize[iHx].XI - (i + 0.5))) / sgg.PML.NumLayers[1][1];
            } else if (i >= SINPML_fullsize[iHx].XE && sgg.PML.NumLayers[1][2] != 0) {
                cm_x[idx] = 1.0 * (i - SINPML_fullsize[iHx].XE + 0.5) / sgg.PML.NumLayers[1][2];
                Icm_x[idx] = 1.0 * (sgg.PML.NumLayers[1][2] - (i - SINPML_fullsize[iHx].XE + 0.5)) / sgg.PML.NumLayers[1][2];
            } else {
                cm_x[idx] = 0.0;
                Icm_x[idx] = 0.0;
            }
        }
        for (int j = hyi; j <= hye; ++j) {
            int idx = j - hyi;
            if (j <= SINPML_fullsize[iHy].YI && sgg.PML.NumLayers[2][1] != 0) {
                ce_y[idx] = 1.0 * (SINPML_fullsize[iHy].YI - j) / sgg.PML.NumLayers[2][1];
                Ice_y[idx] = 1.0 * (sgg.PML.NumLayers[2][1] - (SINPML_fullsize[iHy].YI - j)) / sgg.PML.NumLayers[2][1];
            } else if (j >= SINPML_fullsize[iHy].YE && sgg.PML.NumLayers[2][2] != 0) {
                ce_y[idx] = 1.0 * (j - SINPML_fullsize[iHy].YE) / sgg.PML.NumLayers[2][2];
                Ice_y[idx] = 1.0 * (sgg.PML.NumLayers[2][2] - (j - SINPML_fullsize[iHy].YE)) / sgg.PML.NumLayers[2][2];
            } else {
                ce_y[idx] = 0.0;
                Ice_y[idx] = 0.0;
            }
        }
        for (int j = hyi; j <= hye; ++j) {
            int idx = j - hyi;
            if (j <= SINPML_fullsize[iHy].YI - 1 && sgg.PML.NumLayers[2][1] != 0) {
                cm_y[idx] = 1.0 * (SINPML_fullsize[iHy].YI - (j + 0.5)) / sgg.PML.NumLayers[2][1];
                Icm_y[idx] = 1.0 * (sgg.PML.NumLayers[2][1] - (SINPML_fullsize[iHy].YI - (j + 0.5))) / sgg.PML.NumLayers[2][1];
            } else if (j >= SINPML_fullsize[iHy].YE && sgg.PML.NumLayers[2][2] != 0) {
                cm_y[idx] = 1.0 * (j - SINPML_fullsize[iHy].YE + 0.5) / sgg.PML.NumLayers[2][2];
                Icm_y[idx] = 1.0 * (sgg.PML.NumLayers[2][2] - (j - SINPML_fullsize[iHy].YE + 0.5)) / sgg.PML.NumLayers[2][2];
            } else {
                cm_y[idx] = 0.0;
                Icm_y[idx] = 0.0;
            }
        }
        for (int k = hzi; k <= hze; ++k) {
            int idx = k - hzi;
            if (k <= SINPML_fullsize[iHz].ZI && sgg.PML.NumLayers[3][1] != 0) {
                ce_z[idx] = 1.0 * (SINPML_fullsize[iHz].ZI - k) / sgg.PML.NumLayers[3][1];
                Ice_z[idx] = 1.0 * (sgg.PML.NumLayers[3][1] - (SINPML_fullsize[iHz].ZI - k)) / sgg.PML.NumLayers[3][1];
            } else if (k >= SINPML_fullsize[iHz].ZE && sgg.PML.NumLayers[3][2] != 0) {
                ce_z[idx] = 1.0 * (k - SINPML_fullsize[iHz].ZE) / sgg.PML.NumLayers[3][2];
                Ice_z[idx] = 1.0 * (sgg.PML.NumLayers[3][2] - (k - SINPML_fullsize[iHz].ZE)) / sgg.PML.NumLayers[3][2];
            } else {
                ce_z[idx] = 0.0;
                Ice_z[idx] = 0.0;
            }
        }
        for (int k = hzi; k <= hze; ++k) {
            int idx = k - hzi;
            if (k <= SINPML_fullsize[iHz].ZI - 1 && sgg.PML.NumLayers[3][1] != 0) {
                cm_z[idx] = 1.0 * (SINPML_fullsize[iHz].ZI - (k + 0.5)) / sgg.PML.NumLayers[3][1];
                Icm_z[idx] = 1.0 * (sgg.PML.NumLayers[3][1] - (SINPML_fullsize[iHz].ZI - (k + 0.5))) / sgg.PML.NumLayers[3][1];
            } else if (k >= SINPML_fullsize[iHz].ZE && sgg.PML.NumLayers[3][2] != 0) {
                cm_z[idx] = 1.0 * (k - SINPML_fullsize[iHz].ZE + 0.5) / sgg.PML.NumLayers[3][2];
                Icm_z[idx] = 1.0 * (sgg.PML.NumLayers[3][2] - (k - SINPML_fullsize[iHz].ZE + 0.5)) / sgg.PML.NumLayers[3][2];
            } else {
                cm_z[idx] = 0.0;
                Icm_z[idx] = 0.0;
            }
        }

        calc_cpmlconstants(sgg, Idxe, Idye, Idze, Idxh, Idyh, Idzh, eps0, mu0);

        // Fake coms and ends
        if (!sgg.Border.IsDownPML) {
            for (int i = 1; i <= 6; ++i) PMLc[i].ZI(1) = PMLc[i].ZE(1) + 100;
        }
        if (!sgg.Border.IsUpPML) {
            for (int i = 1; i <= 6; ++i) PMLc[i].ZI(2) = PMLc[i].ZE(2) + 100;
        }
        if (!sgg.Border.IsLeftPML) {
            for (int i = 1; i <= 6; ++i) PMLc[i].ZI(3) = PMLc[i].ZE(3) + 100;
        }
        if (!sgg.Border.IsRightPML) {
            for (int i = 1; i <= 6; ++i) PMLc[i].ZI(4) = PMLc[i].ZE(4) + 100;
        }
        if (!sgg.Border.IsFrontPML) {
            for (int i = 1; i <= 6; ++i) PMLc[i].ZI(6) = PMLc[i].ZE(6) + 100;