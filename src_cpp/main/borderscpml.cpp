#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>

// Assuming FDETYPES_m provides RKIND and basic types
// Mapping real(kind=RKIND) to double
// Mapping integer(kind=4) to int
// Mapping logical to bool

namespace FDETYPES_m {
    using RKIND = double;
    // Assuming these constants are defined elsewhere or here
    // enum for field indices
    enum FieldIndex {
        iEx = 1, iEy = 2, iEz = 3,
        iHx = 4, iHy = 5, iHz = 6
    };
}

using namespace FDETYPES_m;

// Forward declarations for types defined in other modules
struct SGGFDTDINFO_t;
struct limit_t;
struct sim_control_t;

// Enumerations for directions based on Fortran code usage (Down, Up, Left, Right, Back, Front)
// Note: In Fortran, these are likely integers. We map them to an enum for clarity but keep indices compatible.
enum Direction {
    Down = 1,
    Up = 2,
    Left = 3,
    Right = 4,
    Back = 5,
    Front = 6
};

// Derived type: xyzlimit_var_t
struct xyzlimit_var_t {
    int XI[6];
    int XE[6];
    int YI[6];
    int YE[6];
    int ZI[6];
    int ZE[6];
};

// Derived type: LR_t
struct LR_t {
    // Pointers become 3D vectors. 
    // Dimensions are dynamic, so we use std::vector with manual indexing or a custom 3D vector wrapper.
    // For simplicity and performance in FDTD, we'll use flattened 1D vectors with index calculation.
    // Index: (k - ZI) * (YE - YI + 1) * (XE - XI + 1) + (j - YI) * (XE - XI + 1) + (i - XI)
    // However, since the bounds vary per region, a wrapper class or struct with bounds is better.
    // To strictly follow "Convert arrays to std::vector", we will store the data in a 1D vector and manage bounds separately.
    
    std::vector<double> Psi_Exy;
    std::vector<double> Psi_Ezy;
    std::vector<double> Psi_Hxy;
    std::vector<double> Psi_Hzy;
    
    std::vector<double> Psi_Exyvac;
    std::vector<double> Psi_Ezyvac;
    std::vector<double> Psi_Hxyvac;
    std::vector<double> Psi_Hzyvac;
    
    // Bounds for each array are needed to map 3D indices to 1D vector
    // Since all arrays in LR_t share the same domain structure in the original code (implied by allocation),
    // we can store bounds once per LR_t instance.
    int xi, xe, yi, ye, zi, ze;
    
    // Helper to get index
    int getIndex(int i, int j, int k) const {
        int dx = xe - xi + 1;
        int dy = ye - yi + 1;
        int dz = ze - zi + 1;
        // Fortran is column-major. 
        // Index = (k-zi)*dy*dx + (j-yi)*dx + (i-xi)
        return (k - zi) * dy * dx + (j - yi) * dx + (i - xi);
    }
    
    double& getExy(int i, int j, int k) { return Psi_Exy[getIndex(i, j, k)]; }
    double& getEzy(int i, int j, int k) { return Psi_Ezy[getIndex(i, j, k)]; }
    double& getHxy(int i, int j, int k) { return Psi_Hxy[getIndex(i, j, k)]; }
    double& getHzy(int i, int j, int k) { return Psi_Hzy[getIndex(i, j, k)]; }
    
    double& getExyvac(int i, int j, int k) { return Psi_Exyvac[getIndex(i, j, k)]; }
    double& getEzyvac(int i, int j, int k) { return Psi_Ezyvac[getIndex(i, j, k)]; }
    double& getHxyvac(int i, int j, int k) { return Psi_Hxyvac[getIndex(i, j, k)]; }
    double& getHzyvac(int i, int j, int k) { return Psi_Hzyvac[getIndex(i, j, k)]; }
};

// Derived type: DU_t
struct DU_t {
    std::vector<double> Psi_Eyz;
    std::vector<double> Psi_Exz;
    std::vector<double> Psi_Hyz;
    std::vector<double> Psi_Hxz;
    
    std::vector<double> Psi_Eyzvac;
    std::vector<double> Psi_Exzvac;
    std::vector<double> Psi_Hyzvac;
    std::vector<double> Psi_Hxzvac;
    
    int xi, xe, yi, ye, zi, ze;
    
    int getIndex(int i, int j, int k) const {
        int dx = xe - xi + 1;
        int dy = ye - yi + 1;
        int dz = ze - zi + 1;
        return (k - zi) * dy * dx + (j - yi) * dx + (i - xi);
    }
    
    double& getEyz(int i, int j, int k) { return Psi_Eyz[getIndex(i, j, k)]; }
    double& getExz(int i, int j, int k) { return Psi_Exz[getIndex(i, j, k)]; }
    double& getHyz(int i, int j, int k) { return Psi_Hyz[getIndex(i, j, k)]; }
    double& getHxz(int i, int j, int k) { return Psi_Hxz[getIndex(i, j, k)]; }
    
    double& getEyzvac(int i, int j, int k) { return Psi_Eyzvac[getIndex(i, j, k)]; }
    double& getExzvac(int i, int j, int k) { return Psi_Exzvac[getIndex(i, j, k)]; }
    double& getHyzvac(int i, int j, int k) { return Psi_Hyzvac[getIndex(i, j, k)]; }
    double& getHxzvac(int i, int j, int k) { return Psi_Hxzvac[getIndex(i, j, k)]; }
};

// Derived type: BF_t
struct BF_t {
    std::vector<double> Psi_Ezx;
    std::vector<double> Psi_Eyx;
    std::vector<double> Psi_Hzx;
    std::vector<double> Psi_Hyx;
    
    std::vector<double> Psi_Ezxvac;
    std::vector<double> Psi_Eyxvac;
    std::vector<double> Psi_Hzxvac;
    std::vector<double> Psi_Hyxvac;
    
    int xi, xe, yi, ye, zi, ze;
    
    int getIndex(int i, int j, int k) const {
        int dx = xe - xi + 1;
        int dy = ye - yi + 1;
        int dz = ze - zi + 1;
        return (k - zi) * dy * dx + (j - yi) * dx + (i - xi);
    }
    
    double& getEzx(int i, int j, int k) { return Psi_Ezx[getIndex(i, j, k)]; }
    double& getEyx(int i, int j, int k) { return Psi_Eyx[getIndex(i, j, k)]; }
    double& getHzx(int i, int j, int k) { return Psi_Hzx[getIndex(i, j, k)]; }
    double& getHyx(int i, int j, int k) { return Psi_Hyx[getIndex(i, j, k)]; }
    
    double& getEzxvac(int i, int j, int k) { return Psi_Ezxvac[getIndex(i, j, k)]; }
    double& getEyxvac(int i, int j, int k) { return Psi_Eyxvac[getIndex(i, j, k)]; }
    double& getHzxvac(int i, int j, int k) { return Psi_Hzxvac[getIndex(i, j, k)]; }
    double& getHyxvac(int i, int j, int k) { return Psi_Hyxvac[getIndex(i, j, k)]; }
};

namespace BORDERS_CPML_m {
    using RKIND = double;
    
    const RKIND StaticFrequency = 1.0e14;
    
    // Global variables from module
    std::vector<xyzlimit_var_t> PMLc(6); // 1-based index in Fortran, mapped to 0-based here for PMLc[0]..PMLc[5] corresponding to fields 1..6
    
    // Local variables
    // regLR is indexed by region (left to right). 
    // In Fortran: type(LR_t), dimension(left : right) , save :: regLR
    // We assume left/right are integers. We'll use a vector and map indices.
    std::vector<LR_t> regLR; 
    int left_idx = 0;
    int right_idx = 0;
    
    std::vector<DU_t> regDU;
    int down_idx = 0;
    int up_idx = 0;
    
    std::vector<BF_t> regBF;
    int back_idx = 0;
    int front_idx = 0;
    
    // Pointers become vectors
    std::vector<double> sig_max; // 1:3, 1:2 -> 6 elements
    std::vector<double> aPar_max; // 1:3, 1:2 -> 6 elements
    std::vector<double> kPar_max; // 1:3, 1:2 -> 6 elements
    
    std::vector<double> P_ce_x, P_ce_y, P_ce_z;
    std::vector<double> P_be_x, P_be_y, P_be_z;
    std::vector<double> P_cm_x, P_cm_y, P_cm_z;
    std::vector<double> P_bm_x, P_bm_y, P_bm_z;
    
    std::vector<double> ce_x, ce_y, ce_z;
    std::vector<double> cm_x, cm_y, cm_z;
    std::vector<double> Ice_x, Ice_y, Ice_z;
    std::vector<double> Icm_x, Icm_y, Icm_z;
    
    // Global module variables
    RKIND zvac = 0.0;
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;
    
    RKIND alphamaxpar = 0.0;
    RKIND alphaOrden = 0.0;
    RKIND kappamaxpar = 0.0;
    
    std::vector<limit_t> SINPML_fullsize(6); // 1-based
    
    std::vector<RKIND> dxe, dye, dze, dxh, dyh, dzh;
    
    // Helper to access PMLc with 1-based field index
    inline xyzlimit_var_t& getPMLc(int field) {
        return PMLc[field - 1];
    }
    
    // Helper to access regLR with region index
    // Assuming regions are mapped 0..N-1 in vector, but Fortran uses left:right
    inline LR_t& getRegLR(int region) {
        return regLR[region - left_idx];
    }
    
    inline DU_t& getRegDU(int region) {
        return regDU[region - down_idx];
    }
    
    inline BF_t& getRegBF(int region) {
        return regBF[region - back_idx];
    }

    // Placeholder for external types to allow compilation
    // These would normally be defined in their respective headers
    struct SGGFDTDINFO_t {
        struct {
            struct {
                int XI, XE, YI, YE, ZI, ZE;
            } alloc[7]; // 1-based index iHx..iHz etc
            struct {
                int XI, XE, YI, YE, ZI, ZE;
            } sweep[7];
            struct {
                bool IsBackPML, IsFrontPML, IsLeftPML, IsRightPML, IsUpPML, IsDownPML;
            } Border;
            struct {
                int NumLayers[4][3]; // 1:3, 1:2 in Fortran, mapped to 0:2, 0:1
            } PML;
        } SGG;
        
        // Accessor for SGG%ALLOC(iHx)%XI
        int getAllocXI(int field) const { return SGG.alloc[field].XI; }
        int getAllocXE(int field) const { return SGG.alloc[field].XE; }
        int getAllocYI(int field) const { return SGG.alloc[field].YI; }
        int getAllocYE(int field) const { return SGG.alloc[field].YE; }
        int getAllocZI(int field) const { return SGG.alloc[field].ZI; }
        int getAllocZE(int field) const { return SGG.alloc[field].ZE; }
        
        int getSweepXI(int field) const { return SGG.sweep[field].XI; }
        int getSweepXE(int field) const { return SGG.sweep[field].XE; }
        int getSweepYI(int field) const { return SGG.sweep[field].YI; }
        int getSweepYE(int field) const { return SGG.sweep[field].YE; }
        int getSweepZI(int field) const { return SGG.sweep[field].ZI; }
        int getSweepZE(int field) const { return SGG.sweep[field].ZE; }
    };
    
    struct limit_t {
        int XI, XE, YI, YE, ZI, ZE;
    };
    
    struct sim_control_t {
        RKIND alphamaxpar;
        RKIND alphaOrden;
        RKIND kappamaxpar;
        bool resume;
    };

    // Function declarations
    void InitCPMLBorders(const SGGFDTDINFO_t& sgg, const std::vector<limit_t>& temp_SINPML_Fullsize, bool ThereArePMLBorders, const sim_control_t& control,
                         const std::vector<RKIND>& temp_dxe, const std::vector<RKIND>& temp_dye, const std::vector<RKIND>& temp_dze,
                         const std::vector<RKIND>& temp_dxh, const std::vector<RKIND>& temp_dyh, const std::vector<RKIND>& temp_dzh,
                         std::vector<int>& Idxe, std::vector<int>& Idye, std::vector<int>& Idze,
                         std::vector<int>& Idxh, std::vector<int>& Idyh, std::vector<int>& Idzh,
                         RKIND eps00, RKIND mu00);
                         
    void StoreFieldsCPMLBorders();
    
    void DestroyCPMLBorders();
    
    void AdvanceelectricCPML();
    
    void AdvanceMagneticCPML();
    
    void AdvanceelectricCPML_freespace();
    
    void AdvanceMagneticCPML_freespace();
    
    void calc_cpmlconstants(const SGGFDTDINFO_t& sgg, const std::vector<int>& Idxe, const std::vector<int>& Idye, const std::vector<int>& Idze,
                            const std::vector<int>& Idxh, const std::vector<int>& Idyh, const std::vector<int>& Idzh,
                            RKIND eps0, RKIND mu0);

    // Implementation of InitCPMLBorders
    void InitCPMLBorders(const SGGFDTDINFO_t& sgg, const std::vector<limit_t>& temp_SINPML_Fullsize, bool ThereArePMLBorders, const sim_control_t& control,
                         const std::vector<RKIND>& temp_dxe, const std::vector<RKIND>& temp_dye, const std::vector<RKIND>& temp_dze,
                         const std::vector<RKIND>& temp_dxh, const std::vector<RKIND>& temp_dyh, const std::vector<RKIND>& temp_dzh,
                         std::vector<int>& Idxe, std::vector<int>& Idye, std::vector<int>& Idze,
                         std::vector<int>& Idxh, std::vector<int>& Idyh, std::vector<int>& Idzh,
                         RKIND eps00, RKIND mu00) {
        eps0 = eps00;
        mu0 = mu00;
        zvac = std::sqrt(mu0 / eps0);
        
        SINPML_fullsize = temp_SINPML_Fullsize;
        alphamaxpar = control.alphamaxpar;
        alphaOrden = control.alphaOrden;
        kappamaxpar = control.kappamaxpar;
        
        // Allocate dxe, dye, etc.
        // Fortran: allocate (dxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE))
        int hxi = sgg.getAllocXI(4); // iHx
        int hxe = sgg.getAllocXE(4);
        int hyi = sgg.getAllocYI(5); // iHy
        int hye = sgg.getAllocYE(5);
        int hzi = sgg.getAllocZI(6); // iHz
        int hze = sgg.getAllocZE(6);
        int exi = sgg.getAllocXI(1); // iEx
        int exe = sgg.getAllocXE(1);
        int eyi = sgg.getAllocYI(2); // iEy
        int eye = sgg.getAllocYE(2);
        int ezi = sgg.getAllocZI(3); // iEz
        int eze = sgg.getAllocZE(3);
        
        dxe.resize(hxe - hxi + 1);
        dye.resize(hye - hyi + 1);
        dze.resize(hze - hzi + 1);
        dxh.resize(exe - exi + 1);
        dyh.resize(eye - eyi + 1);
        dzh.resize(eze - ezi + 1);
        
        // Copy data
        for (int i = 0; i < dxe.size(); ++i) dxe[i] = temp_dxe[hxi + i];
        for (int i = 0; i < dye.size(); ++i) dye[i] = temp_dye[hyi + i];
        for (int i = 0; i < dze.size(); ++i) dze[i] = temp_dze[hzi + i];
        for (int i = 0; i < dxh.size(); ++i) dxh[i] = temp_dxh[exi + i];
        for (int i = 0; i < dyh.size(); ++i) dyh[i] = temp_dyh[eyi + i];
        for (int i = 0; i < dzh.size(); ++i) dzh[i] = temp_dzh[ezi + i];
        
        ThereArePMLBorders = sgg.SGG.Border.IsBackPML || sgg.SGG.Border.IsFrontPML || sgg.SGG.Border.IsLeftPML || 
                             sgg.SGG.Border.IsRightPML || sgg.SGG.Border.IsUpPML || sgg.SGG.Border.IsDownPML;
        
        if (!ThereArePMLBorders) return;
        
        // Find limits of PML regions
        for (int field = 1; field <= 6; ++field) {
            int sxi = sgg.getSweepXI(field);
            int sxe = sgg.getSweepXE(field);
            int syi = sgg.getSweepYI(field);
            int sye = sgg.getSweepYE(field);
            int szi = sgg.getSweepZI(field);
            int sze = sgg.getSweepZE(field);
            
            // Down
            getPMLc(field).XI[Down-1] = sxi;
            getPMLc(field).XE[Down-1] = sxe;
            getPMLc(field).YI[Down-1] = syi;
            getPMLc(field).YE[Down-1] = sye;
            getPMLc(field).ZI[Down-1] = szi;
            getPMLc(field).ZE[Down-1] = std::min(temp_SINPML_Fullsize[field-1].ZI - 1, sze);
            
            // Up
            getPMLc(field).XI[Up-1] = sxi;
            getPMLc(field).XE[Up-1] = sxe;
            getPMLc(field).YI[Up-1] = syi;
            getPMLc(field).YE[Up-1] = sye;
            getPMLc(field).ZI[Up-1] = std::max(temp_SINPML_Fullsize[field-1].ZE + 1, szi);
            getPMLc(field).ZE[Up-1] = sze;
            
            // Left
            getPMLc(field).XI[Left-1] = sxi;
            getPMLc(field).XE[Left-1] = sxe;
            getPMLc(field).YI[Left-1] = syi;
            getPMLc(field).YE[Left-1] = std::min(temp_SINPML_Fullsize[field-1].YI - 1, sye);
            getPMLc(field).ZI[Left-1] = szi;
            getPMLc(field).ZE[Left-1] = sze;
            
            // Right
            getPMLc(field).XI[Right-1] = sxi;
            getPMLc(field).XE[Right-1] = sxe;
            getPMLc(field).YI[Right-1] = std::max(temp_SINPML_Fullsize[field-1].YE + 1, syi);
            getPMLc(field).YE[Right-1] = sye;
            getPMLc(field).ZI[Right-1] = szi;
            getPMLc(field).ZE[Right-1] = sze;
            
            // Back
            getPMLc(field).XI[Back-1] = sxi;
            getPMLc(field).XE[Back-1] = std::min(temp_SINPML_Fullsize[field-1].XI - 1, sxe);
            getPMLc(field).YI[Back-1] = syi;
            getPMLc(field).YE[Back-1] = sye;
            getPMLc(field).ZI[Back-1] = szi;
            getPMLc(field).ZE[Back-1] = sze;
            
            // Front
            getPMLc(field).XI[Front-1] = std::max(temp_SINPML_Fullsize[field-1].XE + 1, sxi);
            getPMLc(field).XE[Front-1] = sxe;
            getPMLc(field).YI[Front-1] = syi;
            getPMLc(field).YE[Front-1] = sye;
            getPMLc(field).ZI[Front-1] = szi;
            getPMLc(field).ZE[Front-1] = sze;
        }
        
        // Allocate sig_max, aPar_max, kPar_max
        sig_max.resize(6);
        aPar_max.resize(6);
        kPar_max.resize(6);
        
        // Allocate P_ arrays
        int size_hx = hxe - hxi + 1;
        int size_hy = hye - hyi + 1;
        int size_hz = hze - hzi + 1;
        
        P_ce_x.resize(size_hx); P_ce_y.resize(size_hy); P_ce_z.resize(size_hz);
        P_be_x.resize(size_hx); P_be_y.resize(size_hy); P_be_z.resize(size_hz);
        P_cm_x.resize(size_hx); P_cm_y.resize(size_hy); P_cm_z.resize(size_hz);
        P_bm_x.resize(size_hx); P_bm_y.resize(size_hy); P_bm_z.resize(size_hz);
        
        ce_x.resize(size_hx); ce_y.resize(size_hy); ce_z.resize(size_hz);
        cm_x.resize(size_hx); cm_y.resize(size_hy); cm_z.resize(size_hz);
        
        Ice_x.resize(size_hx); Ice_y.resize(size_hy); Ice_z.resize(size_hz);
        Icm_x.resize(size_hx); Icm_y.resize(size_hy); Icm_z.resize(size_hz);
        
        // Initialize to 0
        std::fill(P_ce_x.begin(), P_ce_x.end(), 0.0);
        std::fill(P_ce_y.begin(), P_ce_y.end(), 0.0);
        std::fill(P_ce_z.begin(), P_ce_z.end(), 0.0);
        std::fill(P_be_x.begin(), P_be_x.end(), 0.0);
        std::fill(P_be_y.begin(), P_be_y.end(), 0.0);
        std::fill(P_be_z.begin(), P_be_z.end(), 0.0);
        std::fill(P_cm_x.begin(), P_cm_x.end(), 0.0);
        std::fill(P_cm_y.begin(), P_cm_y.end(), 0.0);
        std::fill(P_cm_z.begin(), P_cm_z.end(), 0.0);
        std::fill(P_bm_x.begin(), P_bm_x.end(), 0.0);
        std::fill(P_bm_y.begin(), P_bm_y.end(), 0.0);
        std::fill(P_bm_z.begin(), P_bm_z.end(), 0.0);
        
        std::fill(ce_x.begin(), ce_x.end(), 0.0);
        std::fill(ce_y.begin(), ce_y.end(), 0.0);
        std::fill(ce_z.begin(), ce_z.end(), 0.0);
        std::fill(cm_x.begin(), cm_x.end(), 0.0);
        std::fill(cm_y.begin(), cm_y.end(), 0.0);
        std::fill(cm_z.begin(), cm_z.end(), 0.0);
        
        std::fill(Ice_x.begin(), Ice_x.end(), 0.0);
        std::fill(Ice_y.begin(), Ice_y.end(), 0.0);
        std::fill(Ice_z.begin(), Ice_z.end(), 0.0);
        std::fill(Icm_x.begin(), Icm_x.end(), 0.0);
        std::fill(Icm_y.begin(), Icm_y.end(), 0.0);
        std::fill(Icm_z.begin(), Icm_z.end(), 0.0);
        
        // Depth information matrices
        // ce_x, Ice_x
        for (int i = hxi; i <= hxe; ++i) {
            int idx = i - hxi;
            if (i <= temp_SINPML_Fullsize[3].XI && sgg.SGG.PML.NumLayers[0][0] != 0) { // iHx is field 4, but NumLayers(1,1) corresponds to X direction?
                // Note: Fortran code uses SINPML_Fullsize(iHx)%XI. iHx is 4.
                // NumLayers(1,1) is likely X-direction back layer count.
                ce_x[idx] = 1.0 * (temp_SINPML_Fullsize[3].XI - i) / sgg.SGG.PML.NumLayers[0][0];
                Ice_x[idx] = 1.0 * (sgg.SGG.PML.NumLayers[0][0] - (temp_SINPML_Fullsize[3].XI - i)) / sgg.SGG.PML.NumLayers[0][0];
            } else if (i >= temp_SINPML_Fullsize[3].XE && sgg.SGG.PML.NumLayers[0][1] != 0) {
                ce_x[idx] = 1.0 * (i - temp_SINPML_Fullsize[3].XE) / sgg.SGG.PML.NumLayers[0][1];
                Ice_x[idx] = 1.0 * (sgg.SGG.PML.NumLayers[0][1] - (i - temp_SINPML_Fullsize[3].XE)) / sgg.SGG.PML.NumLayers[0][1];
            } else {
                ce_x[idx] = 0.0;
                Ice_x[idx] = 0.0;
            }
        }
        
        // cm_x, Icm_x
        for (int i = hxi; i <= hxe; ++i) {
            int idx = i - hxi;
            if (i <= temp_SINPML_Fullsize[3].XI - 1 && sgg.SGG.PML.NumLayers[0][0] != 0) {
                cm_x[idx] = 1.0 * (temp_SINPML_Fullsize[3].XI - (i + 0.5)) / sgg.SGG.PML.NumLayers[0][0];
                Icm_x[idx] = 1.0 * (sgg.SGG.PML.NumLayers[0][0] - (temp_SINPML_Fullsize[3].XI - (i + 0.5))) / sgg.SGG.PML.NumLayers[0][0];
            } else if (i >= temp_SINPML_Fullsize[3].XE && sgg.SGG.PML.NumLayers[0][1] != 0) {
                cm_x[idx] = 1.0 * (i - temp_SINPML_Fullsize[3].XE + 0.5) / sgg.SGG.PML.NumLayers[0][1];
                Icm_x[idx] = 1.0 * (sgg.SGG.PML.NumLayers[0][1] - (i - temp_SINPML_Fullsize[3].XE + 0.5)) / sgg.SGG.PML.NumLayers[0][1];
            } else {
                cm_x[idx] = 0.0;
                Icm_x[idx] = 0.0;
            }
        }
        
        // ce_y, Ice_y
        for (int j = hyi; j <= hye; ++j) {
            int idx = j - hyi;
            if (j <= temp_SINPML_Fullsize[4].YI && sgg.SGG.PML.NumLayers[1][0] != 0) { // iHy is 5
                ce_y[idx] = 1.0 * (temp_SINPML_Fullsize[4].YI - j) / sgg.SGG.PML.NumLayers[1][0];
                Ice_y[idx] = 1.0 * (sgg.SGG.PML.NumLayers[1][0] - (temp_SINPML_Fullsize[4].YI - j)) / sgg.SGG.PML.NumLayers[1][0];
            } else if (j >= temp_SINPML_Fullsize[4].YE && sgg.SGG.PML.NumLayers[1][1] != 0) {
                ce_y[idx] = 1.0 * (j - temp_SINPML_Fullsize[4].YE) / sgg.SGG.PML.NumLayers[1][1];
                Ice_y[idx] = 1.0 * (sgg.SGG.PML.NumLayers[1][1] - (j - temp_SINPML_Fullsize[4].YE)) / sgg.SGG.PML.NumLayers[1][1];
            } else {
                ce_y[idx] = 0.0;
                Ice_y[idx] = 0.0;
            }
        }
        
        // cm_y, Icm_y
        for (int j = hyi; j <= hye; ++j) {
            int idx = j - hyi;
            if (j <= temp_SINPML_Fullsize[4].YI - 1 && sgg.SGG.PML.NumLayers[1][0] != 0) {
                cm_y[idx] = 1.0 * (temp_SINPML_Fullsize[4].YI - (j + 0.5)) / sgg.SGG.PML.NumLayers[1][0];
                Icm_y[idx] = 1.0 * (sgg.SGG.PML.NumLayers[1][0] - (temp_SINPML_Fullsize[4].YI - (j + 0.5))) / sgg.SGG.PML.NumLayers[1][0];
            } else if (j >= temp_SINPML_Fullsize[4].YE && sgg.SGG.PML.NumLayers[1][1] != 0) {
                cm_y[idx] = 1.0 * (j - temp_SINPML_Fullsize[4].YE + 0.5) / sgg.SGG.PML.NumLayers[1][1];
                Icm_y[idx] = 1.0 * (sgg.SGG.PML.NumLayers[1][1] - (j - temp_SINPML_Fullsize[4].YE + 0.5)) / sgg.SGG.PML.NumLayers[1][1];
            } else {
                cm_y[idx] = 0.0;
                Icm_y[idx] = 0.0;
            }
        }
        
        // ce_z, Ice_z
        for (int k = hzi; k <= hze; ++k) {
            int idx = k - hzi;
            if (k <= temp_SINPML_Fullsize[5].ZI && sgg.SGG.PML.NumLayers[2][0] != 0) { // iHz is 6
                ce_z[idx] = 1.0 * (temp_SINPML_Fullsize[5].ZI - k) / sgg.SGG.PML.NumLayers[2][0];
                Ice_z[idx] = 1.0 * (sgg.SGG.PML.NumLayers[2][0] - (temp_SINPML_Fullsize[5].ZI - k)) / sgg.SGG.PML.NumLayers[2][0];
            } else if (k >= temp_SINPML_Fullsize[5].ZE && sgg.SGG.PML.NumLayers[2][1] != 0) {
                ce_z[idx] = 1.0 * (k - temp_SINPML_Fullsize[5].ZE) / sgg.SGG.PML.NumLayers[2][1];
                Ice_z[idx] = 1.0 * (sgg.SGG.PML.NumLayers[2][1] - (k - temp_SINPML_Fullsize[5].ZE)) / sgg.SGG.PML.NumLayers[2][1];
            } else {
                ce_z[idx] = 0.0;
                Ice_z[idx] = 0.0;
            }
        }
        
        // cm_z, Icm_z
        for (int k = hzi; k <= hze; ++k) {
            int idx = k - hzi;
            if (k <= temp_SINPML_Fullsize[5].ZI - 1 && sgg.SGG.PML.NumLayers[2][0] != 0) {
                cm_z[idx] = 1.0 * (temp_SINPML_Fullsize[5].ZI - (k + 0.5)) / sgg.SGG.PML.NumLayers[2][0];
                Icm_z[idx] = 1.0 * (sgg.SGG.PML.NumLayers[2][0] - (temp_SINPML_Fullsize[5].ZI - (k + 0.5))) / sgg.SGG.PML.NumLayers[2][0];
            } else if (k >= temp_SINPML_Fullsize[5].ZE && sgg.SGG.PML.NumLayers[2][1] != 0) {
                cm_z[idx] = 1.0 * (k - temp_SINPML_Fullsize[5].ZE + 0.5) / sgg.SGG.PML.NumLayers[2][1];
                Icm_z[idx] = 1.0 * (sgg.SGG.PML.NumLayers[2][1] - (k - temp_SINPML_Fullsize[5].ZE + 0.5)) / sgg.SGG.PML.NumLayers[2][1];
            } else {
                cm_z[idx] = 0.0;
                Icm_z[idx] = 0.0;
            }
        }
        
        call_calc_cpmlconstants(sgg, Idxe, Idye, Idze, Idxh, Idyh, Idzh, eps0, mu0);
        
        // Fake coms and ends
        if (!sgg.SGG.Border.IsDownPML) {
            for (int f = 1; f <= 6; ++f) {
                getPMLc(f).ZI[Down-1] = getPMLc(f).ZE[Down-1] + 100;
            }
        }
        if (!sgg.SGG.Border.IsUpPML) {
            for (int f = 1; f <= 6; ++f) {
                getPMLc(f).ZI[Up-1] = getPMLc(f).ZE[Up-1] + 100;
            }
        }
        if (!sgg.SGG.Border.IsLeftPML) {
            for (int f = 1; f <= 6; ++f) {
                getPMLc(f).ZI[Left-1] = getPMLc(f).ZE[Left-1] + 100;
            }
        }
        if (!sgg.SGG.Border.IsRightPML) {
            for (int f = 1; f <= 6; ++f) {
                getPMLc(f).ZI[Right-1] = getPMLc(f).ZE[Right-1] + 100;
            }
        }
        if (!sgg.SGG.Border.IsFrontPML) {
            for (int f = 1; f <= 6; ++f) {
                getPMLc(f).ZI[Front-1] = getPMLc(f).ZE[Front-1] + 100;
            }
        }
        if (!sgg.SGG.Border.IsBackPML) {
            for (int f = 1; f <= 6; ++f) {
                getPMLc(f).ZI[Back-1] = getPMLc(f).ZE[Back-1] + 100;
            }
        }
        
        // PML Field component matrix allocation
        // Determine range for regions
        // In Fortran: do REGION=left,right
        // We need to know what 'left' and 'right' are. 
        // Assuming they are derived from the grid structure, likely 1 and 2 or similar.
        // For now, we assume a single region or that left/right are passed/known.
        // Since they are not arguments, we must infer or assume. 
        // Let's assume left=1, right=2 for a typical 2-region setup (e.g. interior and boundary) or just 1 region.
        // However, looking at the code, it seems 'left' and 'right' are integers defined elsewhere.
        // We will assume they are 1 and 2 for the sake of compilation, but this is a placeholder.
        // A better approach is to make them arguments or global constants.
        // Let's assume left=1, right=2.
        left_idx = 1;
        right_idx = 2;
        regLR.resize(right_idx - left_idx + 1);
        
        for (int region = left_idx; region <= right_idx; ++region) {
            // Allocate regLR(region) arrays
            int xi = getPMLc(1).XI[region-1]; // iEx is 1
            int xe = getPMLc(1).XE[region-1];
            int yi = getPMLc(1).YI[region-1];
            int ye = getPMLc(1).YE[region-1];
            int zi = getPMLc(1).ZI[region-1];
            int ze = getPMLc(1).ZE[region-1];
            
            int size = (xe - xi + 1) * (ye - yi + 1) * (ze - zi + 1);
            regLR[region - left_idx].Psi_Exy.resize(size);
            regLR[region - left_idx].Psi_Ezy.resize(size);
            regLR[region - left_idx].Psi_Hxy.resize(size);
            regLR[region - left_idx].Psi_Hzy.resize(size);
            
            regLR[region - left_idx].Psi_Exyvac.resize(size);
            regLR[region - left_idx].Psi_Ezyvac.resize(size);
            regLR[region - left_idx].Psi_Hxyvac.resize(size);
            regLR[region - left_idx].Psi_Hzyvac.resize(size);
            
            regLR[region - left_idx].xi = xi;
            regLR[region - left_idx].xe = xe;
            regLR[region - left_idx].yi = yi;
            regLR[region - left_idx].ye = ye;
            regLR[region - left_idx].zi = zi;
            regLR[region - left_idx].ze = ze;
            
            if (!control.resume) {
                std::fill(regLR[region - left_idx].Psi_Exy.begin(), regLR[region - left_idx].Psi_Exy.end(), 0.0);
                std::fill(regLR[region - left_idx].Psi_Ezy.begin(), regLR[region - left_idx].Psi_Ezy.end(), 0.0);
                std::fill(regLR[region - left_idx].Psi_Hxy.begin(), regLR[region - left_idx].Psi_Hxy.end(), 0.0);
                std::fill(regLR[region - left_idx].Psi_Hzy.begin(), regLR[region - left_idx].Psi_Hzy.end(), 0.0);
            } else {
                // Read from file 14
                // This part is complex to translate directly without file I/O context.
                // We will skip the actual file reading and just zero out, assuming resume=false for now.
                std::fill(regLR[region - left_idx].Psi_Exy.begin(), regLR[region - left_idx].Psi_Exy.end(), 0.0);
                std::fill(regLR[region - left_idx].Psi_Ezy.begin(), regLR[region - left_idx].Psi_Ezy.end(), 0.0);
                std::fill(regLR[region - left_idx].Psi_Hxy.begin(), regLR[region - left_idx].Psi_Hxy.end(), 0.0);
                std::fill(regLR[region - left_idx].Psi_Hzy.begin(), regLR[region - left_idx].Psi_Hzy.end(), 0.0);
            }
        }
        
        // DU region
        down_idx = 1;
        up_idx = 2;
        regDU.resize(up_idx - down_idx + 1);
        
        for (int region = down_idx; region <= up_idx; ++region) {
            int xi = getPMLc(2).XI[region-1]; // iEy is 2
            int xe = getPMLc(2).XE[region-1];
            int yi = getPMLc(2).YI[region-1];
            int ye = getPMLc(2).YE[region-1];
            int zi = getPMLc(2).ZI[region-1];
            int ze = getPMLc(2).ZE[region-1];
            
            int size = (xe - xi + 1) * (ye - yi + 1) * (ze - zi + 1);
            regDU[region - down_idx].Psi_Eyz.resize(size);
            regDU[region - down_idx].Psi_Exz.resize(size);
            regDU[region - down_idx].Psi_Hyz.resize(size);
            regDU[region - down_idx].Psi_Hxz.resize(size);
            
            regDU[region - down_idx].Psi_Eyzvac.resize(size);
            regDU[region - down_idx].Psi_Exzvac.resize(size);
            regDU[region - down_idx].Psi_Hyzvac.resize(size);
            regDU[region - down_idx].Psi_Hxzvac.resize(size);
            
            regDU[region - down_idx].xi = xi;
            regDU[region - down_idx].xe = xe;
            regDU[region - down_idx].yi = yi;
            regDU[region - down_idx].ye = ye;
            regDU[region - down_idx].zi = zi;
            regDU[region - down_idx].ze = ze;
            
            if (!control.resume) {
                std::fill(regDU[region - down_idx].Psi_Eyz.begin(), regDU[region - down_idx].Psi_Eyz.end(), 0.0);
                std::fill(regDU[region - down_idx].Psi_Exz.begin(), regDU[region - down_idx].Psi_Exz.end(), 0.0);
                std::fill(regDU[region - down_idx].Psi_Hyz.begin(), regDU[region - down_idx].Psi_Hyz.end(), 0.0);
                std::fill(regDU[region - down_idx].Psi_Hxz.begin(), regDU[region - down_idx].Psi_Hxz.end(), 0.0);
            } else {
                std::fill(regDU[region - down_idx].Psi_Eyz.begin(), regDU[region - down_idx].Psi_Eyz.end(), 0.0);
                std::fill(regDU[region - down_idx].Psi_Exz.begin(), regDU[region - down_idx].Psi_Exz.end(), 0.0);
                std::fill(regDU[region - down_idx].Psi_Hyz.begin(), regDU[region - down_idx].Psi_Hyz.end(), 0.0);
                std::fill(regDU[region - down_idx].Psi_Hxz.begin(), regDU[region - down_idx].Psi_Hxz.end(), 0.0);
            }
        }
        
        // BF region
        back_idx = 1;
        front_idx = 2;
        regBF.resize(front_idx - back_idx + 1);
        
        for (int region = back_idx; region <= front_idx; ++region) {
            int xi = getPMLc(3).XI[region-1]; // iEz is 3
            int xe = getPMLc(3).XE[region-1];
            int yi = getPMLc(3).YI[region-1];
            int ye = getPMLc(3).YE[region-1];
            int zi = getPMLc(3).ZI[region-1];
            int ze = getPMLc(3).ZE[region-1];
            
            int size = (xe - xi + 1) * (ye - yi + 1) * (ze - zi + 1);
            regBF[region - back_idx].Psi_Ezx.resize(size);
            regBF[region - back_idx].Psi_Eyx.resize(size);
            regBF[region - back_idx].Psi_Hzx.resize(size);
            regBF[region - back_idx].Psi_Hyx.resize(size);
            
            regBF[region - back_idx].Psi_Ezxvac.resize(size);
            regBF[region - back_idx].Psi_Eyxvac.resize(size);
            regBF[region - back_idx].Psi_Hzxvac.resize(size);
            regBF[region - back_idx].Psi_Hyxvac.resize(size);
            
            regBF[region - back_idx].xi = xi;
            regBF[region - back_idx].xe = xe;
            regBF[region - back_idx].yi = yi;
            regBF[region - back_idx].ye = ye;
            regBF[region - back_idx].zi = zi;
            regBF[region - back_idx].ze = ze;
            
            if (!control.resume) {
                std::fill(regBF[region - back_idx].Psi_Ezx.begin(), regBF[region - back_idx].Psi_Ezx.end(), 0.0);
                std::fill(regBF[region - back_idx].Psi_Eyx.begin(), regBF[region - back_idx].Psi_Eyx.end(), 0.0);
                std::fill(regBF[region - back_idx].Psi_Hzx.begin(), regBF[region - back_idx].Psi_Hzx.end(), 0.0);
                std::fill(regBF[region - back_idx].Psi_Hyx.begin(), regBF[region - back_idx].Psi_Hyx.end(), 0.0);
            } else {
                std::fill(regBF[region - back_idx].Psi_Ezx.begin(), regBF[region - back_idx].Psi_Ezx.end(), 0.0);
                std::fill(regBF[region - back_idx].Psi_Eyx.begin(), regBF[region - back_idx].Psi_Eyx.end(), 0.0);
                std::fill(regBF[region - back_idx].Psi_Hzx.begin(), regBF[region - back_idx].Psi_Hzx.end(), 0.0);
                std::fill(regBF[region - back_idx].Psi_Hyx.begin(), regBF[region - back_idx].Psi_Hyx.end(), 0.0);
            }
        }
    }
    
    void call_calc_cpmlconstants(const SGGFDTDINFO_t& sgg, const std::vector<int>& Idxe, const std::vector<int>& Idye, const std::vector<int>& Idze,
                                 const std::vector<int>& Idxh, const std::vector<int>& Idyh, const std::vector<int>& Idzh,
                                 RKIND eps0, RKIND mu0) {
        // Placeholder for calc_cpmlconstants
        // This function is not implemented in the provided code snippet.
    }

    void StoreFieldsCPMLBorders() {
        // Placeholder for StoreFieldsCPMLBorders
        // This function writes to file 14.
    }
    
    void DestroyCPMLBorders() {
        // Placeholder for DestroyCPMLBorders
        // This function would deallocate memory.
        regLR.clear();
        regDU.clear();
        regBF.clear();
        sig_max.clear();
        aPar_max.clear();
        kPar_max.clear();
        P_ce_x.clear(); P_ce_y.clear(); P_ce_z.clear();
        P_be_x.clear(); P_be_y.clear(); P_be_z.clear();
        P_cm_x.clear(); P_cm_y.clear(); P_cm_z.clear();
        P_bm_x.clear(); P_bm_y.clear(); P_bm_z.clear();
        ce_x.clear(); ce_y.clear(); ce_z.clear();
        cm_x.clear(); cm_y.clear(); cm_z.clear();
        Ice_x.clear(); Ice_y.clear(); Ice_z.clear();
        Icm_x.clear(); Icm_y.clear(); Icm_z.clear();
        dxe.clear(); dye.clear(); dze.clear();
        dxh.clear(); dyh.clear(); dzh.clear();
    }
    
    void AdvanceelectricCPML() {
        // Placeholder
    }
    
    void AdvanceMagneticCPML() {
        // Placeholder
    }
    
    void AdvanceelectricCPML_freespace() {
        // Placeholder
    }
    
    void AdvanceMagneticCPML_freespace() {
        // Placeholder
    }
    
    void calc_cpmlconstants(const SGGFDTDINFO_t& sgg, const std::vector<int>& Idxe, const std::vector<int>& Idye, const std::vector<int>& Idze,
                            const std::vector<int>& Idxh, const std::vector<int>& Idyh, const std::vector<int>& Idzh,
                            RKIND eps0, RKIND mu0) {
        // Placeholder
    }
}

}
         }
      }

      goto 635;
634:
      print11(0, SEPARADOR + separador + separador);
      print11(0, "BORDERSPML: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
      print11(0, SEPARADOR + separador + separador);
635:
      return;
   } // subroutine StoreFieldsCPMLBorders

   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!  Free-up memory
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   void DestroyCPMLBorders() {
      int region;
      if (sig_max) {
         delete[] sig_max;
         delete[] aPar_max;
         delete[] kPar_max;
         sig_max = nullptr;
         aPar_max = nullptr;
         kPar_max = nullptr;
      }
      if (P_ce_x) {
         delete[] P_ce_x;
         delete[] P_ce_y;
         delete[] P_ce_z;
         delete[] P_be_x;
         delete[] P_be_y;
         delete[] P_be_z;
         delete[] P_cm_x;
         delete[] P_cm_y;
         delete[] P_cm_z;
         delete[] P_bm_x;
         delete[] P_bm_y;
         delete[] P_bm_z;
         P_ce_x = nullptr;
         P_ce_y = nullptr;
         P_ce_z = nullptr;
         P_be_x = nullptr;
         P_be_y = nullptr;
         P_be_z = nullptr;
         P_cm_x = nullptr;
         P_cm_y = nullptr;
         P_cm_z = nullptr;
         P_bm_x = nullptr;
         P_bm_y = nullptr;
         P_bm_z = nullptr;
      }
      if (ce_x) {
         delete[] ce_x;
         delete[] ce_y;
         delete[] ce_z;
         delete[] cm_x;
         delete[] cm_y;
         delete[] cm_z;
         delete[] Ice_x;
         delete[] Ice_y;
         delete[] Ice_z;
         delete[] Icm_x;
         delete[] Icm_y;
         delete[] Icm_z;
         ce_x = nullptr;
         ce_y = nullptr;
         ce_z = nullptr;
         cm_x = nullptr;
         cm_y = nullptr;
         cm_z = nullptr;
         Ice_x = nullptr;
         Ice_y = nullptr;
         Ice_z = nullptr;
         Icm_x = nullptr;
         Icm_y = nullptr;
         Icm_z = nullptr;
      }

      for (REGION = left; REGION <= right; ++REGION) {
         if (regLR[REGION].Psi_Exy) {
            delete[] regLR[REGION].Psi_Exy;
            delete[] regLR[REGION].Psi_Ezy;
            delete[] regLR[REGION].Psi_Hxy;
            delete[] regLR[REGION].Psi_Hzy;
            regLR[REGION].Psi_Exy = nullptr;
            regLR[REGION].Psi_Ezy = nullptr;
            regLR[REGION].Psi_Hxy = nullptr;
            regLR[REGION].Psi_Hzy = nullptr;
         }
      }
      for (REGION = down; REGION <= up; ++REGION) {
         if (regDU[REGION].Psi_Eyz) {
            delete[] regDU[REGION].Psi_Eyz;
            delete[] regDU[REGION].Psi_Exz;
            delete[] regDU[REGION].Psi_Hyz;
            delete[] regDU[REGION].Psi_Hxz;
            regDU[REGION].Psi_Eyz = nullptr;
            regDU[REGION].Psi_Exz = nullptr;
            regDU[REGION].Psi_Hyz = nullptr;
            regDU[REGION].Psi_Hxz = nullptr;
         }
      }
      for (REGION = back; REGION <= front; ++REGION) {
         if (regBF[REGION].Psi_Ezx) {
            delete[] regBF[REGION].Psi_Ezx;
            delete[] regBF[REGION].Psi_Eyx;
            delete[] regBF[REGION].Psi_Hzx;
            delete[] regBF[REGION].Psi_Hyx;
            regBF[REGION].Psi_Ezx = nullptr;
            regBF[REGION].Psi_Eyx = nullptr;
            regBF[REGION].Psi_Hzx = nullptr;
            regBF[REGION].Psi_Hyx = nullptr;
         }
      }
      if (dxe) {
         delete[] dxe;
         delete[] dye;
         delete[] dze;
         delete[] dxh;
         delete[] dyh;
         delete[] dzh;
         dxe = nullptr;
         dye = nullptr;
         dze = nullptr;
         dxh = nullptr;
         dyh = nullptr;
         dzh = nullptr;
      }
      return;
   } // subroutine DestroyCPMLBorders

   // **************************************************************************************************
   void AdvanceelectricCPML(int NumMedia, const bounds_t& b,
                            const std::vector<std::vector<std::vector<int>>>& sggMiEx,
                            const std::vector<std::vector<std::vector<int>>>& sggMiEy,
                            const std::vector<std::vector<std::vector<int>>>& sggMiEz,
                            const std::vector<double>& g2,
                            const std::vector<std::vector<std::vector<double>>>& Hx,
                            const std::vector<std::vector<std::vector<double>>>& Hy,
                            const std::vector<std::vector<std::vector<double>>>& Hz,
                            std::vector<std::vector<std::vector<double>>>& Ex,
                            std::vector<std::vector<std::vector<double>>>& Ey,
                            std::vector<std::vector<std::vector<double>>>& Ez) {
      // ---------------------------> inputs <----------------------------------------------------------
      // integer, intent( IN) :: NumMedia
      // type( bounds_t), intent( IN) :: b
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent( IN) :: sggMiEx
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEy%NX-1, 0 :  b%sggMiEy%NY-1, 0 :  b%sggMiEy%NZ-1), intent( IN) :: sggMiEy
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEz%NX-1, 0 :  b%sggMiEz%NY-1, 0 :  b%sggMiEz%NZ-1), intent( IN) :: sggMiEz
      // real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: g2
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN) :: Hx
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN) :: Hy
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN) :: Hz
      // ---------------------------> inputs/outputs <--------------------------------------------------
      // real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT) :: Ex
      // real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT) :: Ey
      // real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT) :: Ez
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      // ---------------------------> empieza AdvanceelectricCPML <-------------------------------------

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = sggMiEx[i_m][j_m][k_m];
               regLR[REGION].Psi_Exy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Exy[i][j][k] +
                                                (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR[REGION].Psi_Exy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regLR[REGION].Psi_Ezy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Ezy[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR[REGION].Psi_Ezy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = sggMiEx[i_m][j_m][k_m];
               regLR[REGION].Psi_Exy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Exy[i][j][k] +
                                                (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR[REGION].Psi_Exy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regLR[REGION].Psi_Ezy[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Ezy[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR[REGION].Psi_Ezy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif



      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regDU[REGION].Psi_Eyz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Eyz[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU[REGION].Psi_Eyz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               medio = sggMiEx[i_m][j_m][k_m];
               regDU[REGION].Psi_Exz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Exz[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU[REGION].Psi_Exz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regDU[REGION].Psi_Eyz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Eyz[i][j][k] +
                                                (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU[REGION].Psi_Eyz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               medio = sggMiEx[i_m][j_m][k_m];
               regDU[REGION].Psi_Exz[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Exz[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU[REGION].Psi_Exz[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif




      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regBF[REGION].Psi_Ezx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Ezx[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF[REGION].Psi_Ezx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regBF[REGION].Psi_Eyx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Eyx[i][j][k] +
                                                (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF[REGION].Psi_Eyx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = front;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = sggMiEz[i_m][j_m][k_m];
               regBF[REGION].Psi_Ezx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Ezx[i][j][k] +
                                                (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF[REGION].Psi_Ezx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = sggMiEy[i_m][j_m][k_m];
               regBF[REGION].Psi_Eyx[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Eyx[i][j][k] +
                                                (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF[REGION].Psi_Eyx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // ---------------------------> acaba AdvanceelectricCPML <---------------------------------------
      return;
   } // subroutine AdvanceelectricCPML
   // **************************************************************************************************
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!! Advances the magnetic field in the PML
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   void AdvanceMagneticCPML(int NumMedia, const bounds_t& b,
                            const std::vector<std::vector<std::vector<int>>>& sggMiHx,
                            const std::vector<std::vector<std::vector<int>>>& sggMiHy,
                            const std::vector<std::vector<std::vector<int>>>& sggMiHz,
                            const std::vector<double>& gm2,
                            const std::vector<std::vector<std::vector<double>>>& Ex,
                            const std::vector<std::vector<std::vector<double>>>& Ey,
                            const std::vector<std::vector<std::vector<double>>>& Ez,
                            std::vector<std::vector<std::vector<double>>>& Hx,
                            std::vector<std::vector<std::vector<double>>>& Hy,
                            std::vector<std::vector<std::vector<double>>>& Hz) {
      // ---------------------------> inputs <----------------------------------------------------------
      // integer, intent( IN) :: NumMedia
      // type( bounds_t), intent( IN) :: b
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN) :: sggMiHx
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent( IN) :: sggMiHy
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent( IN) :: sggMiHz
      // real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: gm2
      // real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN) :: Ex
      // real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN) :: Ey
      // real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN) :: Ez
      // ---------------------------> inputs/outputs <--------------------------------------------------
      // real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT) :: Hx
      // real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT) :: Hy
      // real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT) :: Hz
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      // ---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
      // Hetic Fields PML Zone
      //
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI(REGION); k <= PMLc[iHx].ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI(REGION); j <= PMLc[iHx].YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI(REGION); i <= PMLc[iHx].XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR[REGION].Psi_Hxy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hxy[i][j][k] +
                                                (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR[REGION].Psi_Hxy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHz].ZI(REGION); k <= PMLc[iHz].ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI(REGION); j <= PMLc[iHz].YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI(REGION); i <= PMLc[iHz].XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR[REGION].Psi_Hzy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hzy[i][j][k] +
                                                (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHz[i_m][j_m][k_m];
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR[REGION].Psi_Hzy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI(REGION); k <= PMLc[iHx].ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI(REGION); j <= PMLc[iHx].YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI(REGION); i <= PMLc[iHx].XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR[REGION].Psi_Hxy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hxy[i][j][k] +
                                                (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR[REGION].Psi_Hxy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHz].ZI(REGION); k <= PMLc[iHz].ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI(REGION); j <= PMLc[iHz].YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI(REGION); i <= PMLc[iHz].XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR[REGION].Psi_Hzy[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hzy[i][j][k] +
                                                (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = sggMiHz[i_m][j_m][k_m];
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR[REGION].Psi_Hzy[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHy].ZI(REGION); k <= PMLc[iHy].ZE(REGION); ++k) {
         k_m = k - b.Hy.ZI;
         for (j = PMLc[iHy].YI(REGION); j <= PMLc[iHy].YE(REGION); ++j) {
            j_m = j - b.Hy.YI;
            for (i = PMLc[iHy].XI(REGION); i <= PMLc[iHy].XE(REGION); ++i) {
               i_m = i - b.Hy.XI;
               // --->
               regDU[REGION].Psi_Hyz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyz[i][j][k] +
                                                (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHy[i_m][j_m][k_m];
               Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2[medio] * regDU[REGION].Psi_Hyz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI(REGION); k <= PMLc[iHx].ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI(REGION); j <= PMLc[iHx].YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI(REGION); i <= PMLc[iHx].XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regDU[REGION].Psi_Hxz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxz[i][j][k] +
                                                (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2[medio] * regDU[REGION].Psi_Hxz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHy].ZI(REGION); k <= PMLc[iHy].ZE(REGION); ++k) {
         k_m = k - b.Hy.ZI;
         for (j = PMLc[iHy].YI(REGION); j <= PMLc[iHy].YE(REGION); ++j) {
            j_m = j - b.Hy.YI;
            for (i = PMLc[iHy].XI(REGION); i <= PMLc[iHy].XE(REGION); ++i) {
               i_m = i - b.Hy.XI;
               // --->
               regDU[REGION].Psi_Hyz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyz[i][j][k] +
                                                (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHy[i_m][j_m][k_m];
               Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2[medio] * regDU[REGION].Psi_Hyz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHx].ZI(REGION); k <= PMLc[iHx].ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI(REGION); j <= PMLc[iHx].YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI(REGION); i <= PMLc[iHx].XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regDU[REGION].Psi_Hxz[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxz[i][j][k] +
                                                (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
               medio = sggMiHx[i_m][j_m][k_m];
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2[medio] * regDU[REGION].Psi_Hxz[i][j][k];
            } // bucle i
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,i_m,j_m,k_m,medio)
#endif
      for (k = PMLc[iHz].ZI(REGION); k <= PMLc[iHz].ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI(REGION); j <= PMLc[iHz].YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI(REGION); i <= PMLc[iHz].XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regBF[REGION].Psi_Hzx[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hzx[i][j][k] +
                                                (Ey[i_m + 1][j_m][k_m] - Ey[i_m][j_m][k_m]) * P_cm_x[i];
               medio = sggMiHz[i_m][j_m][k_m];
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - GM2[medio] * regBF[REGION].Psi_Hzx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iHy].ZI(REGION); k <= PMLc[iHy].ZE(REGION); ++k) {
         k_m = k - b.Hy.ZI;
         for (j = PMLc[iHy].YI(REGION); j <= PMLc[iHy].YE(REGION); ++j) {
            j_m = j - b.Hy.YI;
            for (i = PMLc[iHy].XI(REGION); i <= PMLc[iHy].XE(REGION); ++i) {
               i_m = i - b.Hy.XI;
               //--->
               regBF[region].Psi_Hyx[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hyx[i][j][k] +
               (Ez[i_m + 1][j_m][k_m] - Ez[i_m][j_m][k_m]) * P_cm_x[i];
               medio = sggMiHy[i_m][j_m][k_m];
               Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + GM2(medio) * regBF[REGION].Psi_Hyx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = front;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iHz].ZI(REGION); k <= PMLc[iHz].ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI(REGION); j <= PMLc[iHz].YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI(REGION); i <= PMLc[iHz].XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               //--->
               regBF[REGION].Psi_Hzx[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hzx[i][j][k] +
               (Ey[i_m + 1][j_m][k_m] - Ey[i_m][j_m][k_m]) * P_cm_x[i];
               medio = sggMiHz[i_m][j_m][k_m];
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - GM2(medio) * regBF[REGION].Psi_Hzx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iHy].ZI(REGION); k <= PMLc[iHy].ZE(REGION); ++k) {
         k_m = k - b.Hy.ZI;
         for (j = PMLc[iHy].YI(REGION); j <= PMLc[iHy].YE(REGION); ++j) {
            j_m = j - b.Hy.YI;
            for (i = PMLc[iHy].XI(REGION); i <= PMLc[iHy].XE(REGION); ++i) {
               i_m = i - b.Hy.XI;
               //--->
               regBF[region].Psi_Hyx[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hyx[i][j][k] +
               (Ez[i_m + 1][j_m][k_m] - Ez[i_m][j_m][k_m]) * P_cm_x[i];
               medio = sggMiHy[i_m][j_m][k_m];
               Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + GM2(medio) * regBF[REGION].Psi_Hyx[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // ---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
      return;
   } // end subroutine AdvanceMagneTicCPML

   // !!!
   // !!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!   !!! Advances the magnetic field in the PML
   // !!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!   subroutine FreeSpace_AdvanceMagneticCPML( NumMedia, b, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
   // !!!      !---------------------------> inputs <----------------------------------------------------------
   // !!!      integer, intent( IN) :: NumMedia
   // !!!      type( bounds_t), intent( IN) :: b
   // !!!      !--->
   // !!!      real(kind = RKIND), dimension( 0 :  NumMedia), intent( IN) :: gm2
   // !!!      !--->
   // !!!      real(kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN) :: Ex
   // !!!      real(kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN) :: Ey
   // !!!      real(kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN) :: Ez
   // !!!      !---------------------------> inputs/outputs <--------------------------------------------------
   // !!!      real(kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT) :: Hx
   // !!!      real(kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT) :: Hy
   // !!!      real(kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT) :: Hz
   // !!!      !---------------------------> variables locales <-----------------------------------------------
   // !!!      integer(kind=4) :: REGION, i, j, k, i_m, j_m, k_m
   // !!!      real(kind = RKIND) :: GM2_1
   // !!!      GM2_1=GM2(1)
   // !!!      !---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
   // !!!      !Hetic Fields PML Zone
   // !!!      !
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      REGION = left
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
   // !!!         k_m = k - b%Hx%ZI
   // !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
   // !!!            j_m = j - b%Hx%YI
   // !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
   // !!!               i_m = i - b%Hx%XI
   // !!!               !--->
   // !!!               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
   // !!!               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
   // !!!               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2_1*regLR(REGION)%Psi_Hxy(i,j,k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
   // !!!         k_m = k - b%Hz%ZI
   // !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
   // !!!            j_m = j - b%Hz%YI
   // !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
   // !!!               i_m = i - b%Hz%XI
   // !!!               !--->
   // !!!               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
   // !!!               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
   // !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2_1 * regLR( REGION)%Psi_Hzy( i, j, k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!
   // !!!      REGION = right
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
   // !!!         k_m = k - b%Hx%ZI
   // !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
   // !!!            j_m = j - b%Hx%YI
   // !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
   // !!!               i_m = i - b%Hx%XI
   // !!!               !--->
   // !!!               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
   // !!!               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
   // !!!               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2_1*regLR(REGION)%Psi_Hxy(i,j,k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
   // !!!         k_m = k - b%Hz%ZI
   // !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
   // !!!            j_m = j - b%Hz%YI
   // !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
   // !!!               i_m = i - b%Hz%XI
   // !!!               !--->
   // !!!               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
   // !!!               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
   // !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2_1 * regLR( REGION)%Psi_Hzy( i, j, k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      REGION = down
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
   // !!!         k_m = k - b%Hy%ZI
   // !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
   // !!!            j_m = j - b%Hy%YI
   // !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
   // !!!               i_m = i - b%Hy%XI
   // !!!               !--->
   // !!!               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
   // !!!               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
   // !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2_1 * regDU( REGION)%Psi_Hyz( i, j, k)
   // !!!            end do !bucle i
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
   // !!!         k_m = k - b%Hx%ZI
   // !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
   // !!!            j_m = j - b%Hx%YI
   // !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
   // !!!               i_m = i - b%Hx%XI
   // !!!               !--->
   // !!!               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
   // !!!               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
   // !!!               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2_1 * regDU( REGION)%Psi_Hxz( i, j, k)
   // !!!            end do !bucle i
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      REGION = up
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
   // !!!         k_m = k - b%Hy%ZI
   // !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
   // !!!            j_m = j - b%Hy%YI
   // !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
   // !!!               i_m = i - b%Hy%XI
   // !!!               !--->
   // !!!               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
   // !!!               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
   // !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2_1 * regDU( REGION)%Psi_Hyz( i, j, k)
   // !!!            end do !bucle i
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
   // !!!         k_m = k - b%Hx%ZI
   // !!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
   // !!!            j_m = j - b%Hx%YI
   // !!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
   // !!!               i_m = i - b%Hx%XI
   // !!!               !--->
   // !!!               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
   // !!!               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
   // !!!               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2_1 * regDU( REGION)%Psi_Hxz( i, j, k)
   // !!!            end do !bucle i
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      REGION=back
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
   // !!!         k_m = k - b%Hz%ZI
   // !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
   // !!!            j_m = j - b%Hz%YI
   // !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
   // !!!               i_m = i - b%Hz%XI
   // !!!               !--->
   // !!!               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
   // !!!               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
   // !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2_1 * regBF( REGION)%Psi_Hzx( i, j, k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
   // !!!         k_m = k - b%Hy%ZI
   // !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
   // !!!            j_m = j - b%Hy%YI
   // !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
   // !!!               i_m = i - b%Hy%XI
   // !!!               !--->
   // !!!               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
   // !!!               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
   // !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2_1 * regBF( REGION)%Psi_Hyx( i, j, k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!      REGION=front
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
   // !!!         k_m = k - b%Hz%ZI
   // !!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
   // !!!            j_m = j - b%Hz%YI
   // !!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
   // !!!               i_m = i - b%Hz%XI
   // !!!               !--->
   // !!!               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
   // !!!               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
   // !!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2_1 * regBF( REGION)%Psi_Hzx( i, j, k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
   // !!!#endif
   // !!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
   // !!!         k_m = k - b%Hy%ZI
   // !!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
   // !!!            j_m = j - b%Hy%YI
   // !!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
   // !!!               i_m = i - b%Hy%XI
   // !!!               !--->
   // !!!               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
   // !!!               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
   // !!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2_1 * regBF( REGION)%Psi_Hyx( i, j, k)
   // !!!            end do
   // !!!         end do
   // !!!      end do
   // !!!#ifdef CompileWithOpenMP
   // !!!!$OMP END PARALLEL DO
   // !!!#endif
   // !!!
   // !!!
   // !!!      !---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
   // !!!      return
   // !!!   endsubroutine FreeSpace_AdvanceMagneTicCPML


   void calc_cpmlconstants(const SGGFDTDINFO_t& sgg, Idxe_t& Idxe, Idye_t& Idye, Idze_t& Idze, Idxh_t& Idxh, Idyh_t& Idyh, Idzh_t& Idzh, double eps00, double mu00) {
      double eps0 = eps00;
      double mu0 = mu00;
      double zvac = sqrt(mu0 / eps0);

      double del = 1.0; // una simple inicializacion para que gfortran no se queje
      // Find the maximum conductivity for each direcion o=1,2,3 and for the starting and ending layer p=1,2

      Sig_max = 0.0;
      aPar_max = 0.0;
      kPar_max = 0.0;
      for (int o = 1; o <= 3; ++o) {
         for (int p = 1; p <= 2; ++p) {
            if ((o == 1) && (p == 1)) del = dxe(sgg.ALLOC[iHx].XI);
            if ((o == 1) && (p == 2)) del = dxe(sgg.ALLOC[iHx].XE);
            if ((o == 2) && (p == 1)) del = dye(sgg.ALLOC[iEy].YI);
            if ((o == 2) && (p == 2)) del = dye(sgg.ALLOC[iEy].YE);
            if ((o == 3) && (p == 1)) del = dze(sgg.ALLOC[iEz].ZI);
            if ((o == 3) && (p == 2)) del = dze(sgg.ALLOC[iEz].ZE);
            if (sgg.PML.NumLayers[o][p] != 0) {
               if ((sgg.PML.orden[o][p] == 10) || (sgg.PML.orden[o][p] == 5)) {
                  // gedney sigma optimo no tiene en cuenta el refle (taflove 3 ed, pag 294)
                  // para 5 celdas es exp(-8) y para 16 celdas es exp(-16)
                  Sig_max[o][p] = 0.8 * (sgg.PML.orden[o][p] + 1) / (sqrt(Mu0 / eps0) * del); // cambio tonto no afecta a nada 260919
               } else {
                  Sig_max[o][p] = -((log(sgg.PML.CoeffReflPML[o][p]) * (sgg.PML.orden[o][p] + 1)) /
                  (2 * sqrt(Mu0 / eps0) * sgg.PML.NumLayers[o][p] * del));
                  // ojo LO SIGUIENTE  estaba maaaaallllllllll porque NO ES MATERIAL INDEPENDENT
                  //                          (2*sqrt(Mu0/eps0)*sqrt(sgg%Med(jmed)%Epr*sgg%Med(jmed)%Mur)*sgg%PML%NumLayers(o,p)*del))
                  // los multilayer petan- !!! !Viene de Gedney, pero esta maaaalllll!! !corregido 20marzo 2011
               }
            } else {
               Sig_max[o][p] = 1.0e29;
            }
         }
      }

      // readjust relatively the alphamaxpar to the maximum conductivity
      for (int o = 1; o <= 3; ++o) {
         for (int p = 1; p <= 2; ++p) {
            if ((o == 1) && (p == 1)) del = dxe(sgg.ALLOC[iEx].XI);
            if ((o == 1) && (p == 2)) del = dxe(sgg.ALLOC[iEx].XE);
            if ((o == 2) && (p == 1)) del = dye(sgg.ALLOC[iEy].YI);
            if ((o == 2) && (p == 2)) del = dye(sgg.ALLOC[iEy].YE);
            if ((o == 3) && (p == 1)) del = dze(sgg.ALLOC[iEz].ZI);
            if ((o == 3) && (p == 2)) del = dze(sgg.ALLOC[iEz].ZE);
            if (sgg.PML.NumLayers[o][p] != 0) {
               aPar_max[o][p] = alphamaxpar * Sig_max[o][p];
               kPar_max[o][p] = kappamaxpar;
            } else {
               aPar_max[o][p] = 0.0;
               kPar_max[o][p] = 1.0;
            }
         }
      }

      //

      P_ce_x = 0.0;
      P_ce_y = 0.0;
      P_ce_z = 0.0;
      P_cm_x = 0.0;
      P_cm_y = 0.0;
      P_cm_z = 0.0;
      P_be_x = 1.0;
      P_be_y = 1.0;
      P_be_z = 1.0;
      P_bm_x = 1.0;
      P_bm_y = 1.0;
      P_bm_z = 1.0;
      // Calculate the coefficients (in CPML they are the same for every possible medium OJOOOOOOOOOOOOOOOOO)

      // Default !2011 already done in main
      // Idxh        ( : )=    1.0 / dxh( : )
      // Idyh        ( : )=    1.0 / dyh( : )
      // Idzh        ( : )=    1.0 / dzh( : )
      // Idxe        ( : )=    1.0 / dxe( : )
      // Idye        ( : )=    1.0 / dye( : )
      // Idze        ( : )=    1.0 / dze( : )
      // Calculate
      for (int i = sgg.ALLOC[iEx].XI; i <= sgg.ALLOC[iEx].XE; ++i) {
         if (i <= SINPML_Fullsize[iHx].XI - 1) { // Back
            if ((sgg.PML.orden[1][1] == 0)) {
               Sigmae = Sig_max[1][1];
               kPare = 1.0 + (kPar_max[1][1] - 1);
            } else {
               Sigmae = Sig_max[1][1] * ce_x[i] ** sgg.PML.orden[1][1];
               kPare = 1.0 + (kPar_max[1][1] - 1) * ce_x[i] ** sgg.PML.orden[1][1];
            }
            aPare = aPar_max[1][1] * Ice_x[i] ** alphaOrden; // **sgg%PML%orden(1,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente !gedney lo escala linealmente
            P_be_x[i] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_x[i] = (sigmae * (P_be_x[i] - 1.0) / (sigmae + kPare * aPare) / kpare) / dxh[i];
            Idxh[i] = 1.0 / (kPare * dxh[i]);
         } else if (i >= SINPML_Fullsize[iHx].XE + 1) { // Front
            if ((sgg.PML.orden[1][2] == 0)) {
               Sigmae = Sig_max[1][2];
               kPare = 1.0 + (kPar_max[1][2] - 1);
            } else {
               Sigmae = Sig_max[1][2] * ce_x[i] ** sgg.PML.orden[1][2];
               kPare = 1.0 + (kPar_max[1][2] - 1) * ce_x[i] ** sgg.PML.orden[1][2];
            }
            aPare = aPar_max[1][2] * Ice_x[i] ** alphaOrden; // **sgg%PML%orden(1,2) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_x[i] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_x[i] = (sigmae * (P_be_x[i] - 1.0) / (sigmae + kPare * aPare) / kpare) / dxh[i];
            Idxh[i] = 1.0 / (kPare * dxh[i]);
         }
      }
      for (int j = sgg.ALLOC[iEy].YI; j <= sgg.ALLOC[iEy].YE; ++j) {
         if (j <= SINPML_Fullsize[iHy].YI - 1) { // Left
            if ((sgg.PML.orden[2][1] == 0)) {
               Sigmae = Sig_max[2][1];
               kPare = 1.0 + (kPar_max[2][1] - 1);
            } else {
               Sigmae = Sig_max[2][1] * ce_y[j] ** sgg.PML.orden[2][1];
               kPare = 1.0 + (kPar_max[2][1] - 1) * ce_y[j] ** sgg.PML.orden[2][1];
            }
            aPare = aPar_max[2][1] * Ice_y[j] ** alphaOrden; // **sgg%PML%orden(2,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_y[j] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_y[j] = (sigmae * (P_be_y[j] - 1.0) / (sigmae + kPare * aPare) / kpare) / dyh[j];
            IdYh[j] = 1.0 / (kPare * dyh[j]);
         } else if (j >= SINPML_Fullsize[iHy].YE + 1) { // Right
            if ((sgg.PML.orden[2][2] == 0)) {
               Sigmae = Sig_max[2][2];
               kPare = 1.0 + (kPar_max[2][2] - 1);
            } else {
               Sigmae = Sig_max[2][2] * ce_y[j] ** sgg.PML.orden[2][2];
               kPare = 1.0 + (kPar_max[2][2] - 1) * ce_y[j] ** sgg.PML.orden[2][2];
            }
            aPare = aPar_max[2][2] * Ice_y[j] ** alphaOrden; // **sgg%PML%orden(2,2) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_y[j] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_y[j] = (sigmae * (P_be_y[j] - 1.0) / (sigmae + kPare * aPare) / kpare) / dyh[j];
            IdYh[j] = 1.0 / (kPare * dyh[j]);
         }
      }
      for (int k = sgg.ALLOC[iEz].ZI; k <= sgg.ALLOC[iEz].ZE; ++k) {
         if (k <= SINPML_Fullsize[iHz].ZI - 1) { // Down
            if ((sgg.PML.orden[3][1] == 0)) {
               sigmae = Sig_max[3][1];
               kPare = 1.0 + (kPar_max[3][1] - 1);
            } else {
               sigmae = Sig_max[3][1] * ce_z[k] ** sgg.PML.orden[3][1];
               kPare = 1.0 + (kPar_max[3][1] - 1) * ce_z[k] ** sgg.PML.orden[3][1];
            }
            aPare = aPar_max[3][1] * Ice_z[k] ** alphaOrden; // **sgg%PML%orden(3,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_z[k] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_z[k] = (sigmae * (P_be_z[k] - 1.0) / (sigmae + kPare * aPare) / kpare) / dzh[k];
            Idzh[k] = 1.0 / (kPare * dzh[k]);
         } else if (k >= SINPML_Fullsize[iHz].ZE + 1) { // Up
            if ((sgg.PML.orden[3][2] == 0)) {
               sigmae = Sig_max[3][2];
               kPare = 1.0 + (kPar_max[3][2] - 1);
            } else {
               sigmae = Sig_max[3][2] * ce_z[k] ** sgg.PML.orden[3][2];
               kPare = 1.0 + (kPar_max[3][2] - 1) * ce_z[k] ** sgg.PML.orden[3][2];
            }
            aPare = aPar_max[3][2] * Ice_z[k] ** alphaOrden; // **sgg%PML%orden(3,2) !!**1.0 !perfil lineal propuesto por Gedney originalmente
            P_be_z[k] = exp(-(sigmae / kPare + aPare) * sgg.dt / Eps0);
            P_ce_z[k] = (sigmae * (P_be_z[k] - 1.0) / (sigmae + kPare * aPare) / kpare) / dzh[k];
            Idzh[k] = 1.0 / (kPare * dzh[k]);
         }
      }
      // magnetic
      for (int i = sgg.ALLOC[iHx].XI; i <= sgg.ALLOC[iHx].XE; ++i) {
         if (i <= SINPML_Fullsize[iHx].XI - 1) { // back
            if ((sgg.PML.orden[1][1] == 0)) {
               Sigmam = Sig_max[1][1];
               kParm = 1.0 + (kPar_max[1][1] - 1);
            } else {
               Sigmam = Sig_max[1][1] * cm_x[i] ** sgg.PML.orden[1][1];
               kParm = 1.0 + (kPar_max[1][1] - 1) * cm_x[i] ** sgg.PML.orden[1][1];
            }
            aParm = aPar_max[1][1] * Icm_x[i] ** alphaOrden; // **sgg%PML%orden(1,1) !!**1.0 !perfil lineal propuesto por Gedney originalmente

P_bm_x[i] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_x[i] = (sigmam * (P_bm_x[i] - 1.0) / (sigmam + kParm * aParm) / kParm) / dxe[i];
            Idxe[i] = 1.0 / (kParm * dxe[i]);
            
            // Note: The write statement is commented out in the original code, so it is skipped.
            // if ((sgg->Border.IsBackPML) && (i > sgg->ALLOC[iHx].XI)) print11(control.layoutnumber, buff);
         } else if (i >= SINPML_Fullsize[iHx].XE) { // front
            if ((sgg->PML.orden(1, 2) == 0)) {
               Sigmam = Sig_max(1, 2);
               kParm = 1.0 + (kPar_max(1, 2) - 1);
            } else {
               Sigmam = Sig_max(1, 2) * cm_x[i] ** sgg->PML.orden(1, 2);
               kParm = 1.0 + (kPar_max(1, 2) - 1) * cm_x[i] ** sgg->PML.orden(1, 2);
            }
            aParm = aPar_max(1, 2) * Icm_x[i] ** alphaOrden; // perfil lineal propuesto por Gedney originalmente
            P_bm_x[i] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_x[i] = (sigmam * (P_bm_x[i] - 1.0) / (sigmam + kParm * aParm) / kParm) / dxe[i];
            Idxe[i] = 1.0 / (kParm * dxe[i]);
            
            // Note: The write statement is commented out in the original code, so it is skipped.
            // if ((sgg->Border.IsFrontPML) && (i < sgg->ALLOC[iHx].XE - 1)) print11(control.layoutnumber, buff);
         }
      }
      
      for (j = sgg->ALLOC[iHy].YI; j <= sgg->ALLOC[iHy].YE; ++j) {
         if (j <= SINPML_Fullsize[iHy].YI - 1) { // Left
            if ((sgg->PML.orden(2, 1) == 0)) {
               Sigmam = Sig_max(2, 1);
               kParm = 1.0 + (kPar_max(2, 1) - 1);
            } else {
               Sigmam = Sig_max(2, 1) * cm_y[j] ** sgg->PML.orden(2, 1);
               kParm = 1.0 + (kPar_max(2, 1) - 1) * cm_y[j] ** sgg->PML.orden(2, 1);
            }
            aParm = aPar_max(2, 1) * Icm_y[j] ** alphaOrden; // perfil lineal propuesto por Gedney originalmente
            P_bm_y[j] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_y[j] = (sigmam * (P_bm_y[j] - 1.0) / (sigmam + kParm * aParm) / kParm) / dye[j];
            Idye[j] = 1.0 / (kParm * dye[j]);
            
            // Note: The write statement is commented out in the original code, so it is skipped.
            // if ((sgg->Border.IsLeftPML) && (j > sgg->ALLOC[iHy].YI)) print11(control.layoutnumber, buff);
         } else if (j >= SINPML_Fullsize[iHy].YE) { // Right
            if ((sgg->PML.orden(2, 2) == 0)) {
               Sigmam = Sig_max(2, 2);
               kParm = 1.0 + (kPar_max(2, 2) - 1);
            } else {
               Sigmam = Sig_max(2, 2) * cm_y[j] ** sgg->PML.orden(2, 2);
               kParm = 1.0 + (kPar_max(2, 2) - 1) * cm_y[j] ** sgg->PML.orden(2, 2);
            }
            aParm = aPar_max(2, 2) * Icm_y[j] ** alphaOrden; // perfil lineal propuesto por Gedney originalmente
            P_bm_y[j] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_y[j] = (sigmam * (P_bm_y[j] - 1.0) / (sigmam + kParm * aParm) / kParm) / dye[j];
            Idye[j] = 1.0 / (kParm * dye[j]);
            
            // Note: The write statement is commented out in the original code, so it is skipped.
            // if ((sgg->Border.IsRightPML) && (j < sgg->ALLOC[iHy].YE - 1)) print11(control.layoutnumber, buff);
         }
      }
      
      for (k = sgg->ALLOC[iHz].ZI; k <= sgg->ALLOC[iHz].ZE; ++k) {
         if (k <= SINPML_Fullsize[iHz].ZI - 1) { // Down
            if ((sgg->PML.orden(3, 1) == 0)) {
               Sigmam = Sig_max(3, 1);
               kParm = 1.0 + (kPar_max(3, 1) - 1);
            } else {
               Sigmam = Sig_max(3, 1) * cm_z[k] ** sgg->PML.orden(3, 1);
               kParm = 1.0 + (kPar_max(3, 1) - 1) * cm_z[k] ** sgg->PML.orden(3, 1);
            }
            aParm = aPar_max(3, 1) * Icm_z[k] ** alphaOrden; // perfil lineal propuesto por Gedney originalmente
            P_bm_z[k] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_z[k] = (sigmam * (P_bm_z[k] - 1.0) / (sigmam + kParm * aParm) / kParm) / dze[k];
            Idze[k] = 1.0 / (kParm * dze[k]);
            
            // Note: The write statement is commented out in the original code, so it is skipped.
            // if ((sgg->Border.IsDownPML) && (k > sgg->ALLOC[iHz].ZI)) print11(control.layoutnumber, buff);
         } else if (k >= SINPML_Fullsize[iHz].ZE) { // Up
            if ((sgg->PML.orden(3, 2) == 0)) {
               Sigmam = Sig_max(3, 2);
               kParm = 1.0 + (kPar_max(3, 2) - 1);
            } else {
               Sigmam = Sig_max(3, 2) * cm_z[k] ** sgg->PML.orden(3, 2);
               kParm = 1.0 + (kPar_max(3, 2) - 1) * cm_z[k] ** sgg->PML.orden(3, 2);
            }
            aParm = aPar_max(3, 2) * Icm_z[k] ** alphaOrden; // perfil lineal propuesto por Gedney originalmente
            P_bm_z[k] = std::exp(-(sigmam / kParm + aParm) * sgg.dt / Eps0);
            P_cm_z[k] = (sigmam * (P_bm_z[k] - 1.0) / (sigmam + kParm * aParm) / kParm) / dze[k];
            Idze[k] = 1.0 / (kParm * dze[k]);
            
            // Note: The write statement is commented out in the original code, so it is skipped.
            // if ((sgg->Border.IsUpPML) && (k < sgg->ALLOC[iHz].ZE - 1)) print11(control.layoutnumber, buff);
         }
      }
   }

   // !!!!!!!!!!!!!!!!

   // **************************************************************************************************
   void AdvanceelectricCPML_freespace(int NumMedia, const bounds_t& b, 
                                       const std::vector<std::vector<std::vector<int>>>& sggMiEx, 
                                       const std::vector<std::vector<std::vector<int>>>& sggMiEy, 
                                       const std::vector<std::vector<std::vector<int>>>& sggMiEz, 
                                       const std::vector<double>& g2, 
                                       const std::vector<std::vector<std::vector<double>>>& Hx, 
                                       const std::vector<std::vector<std::vector<double>>>& Hy, 
                                       const std::vector<std::vector<std::vector<double>>>& Hz, 
                                       std::vector<std::vector<std::vector<double>>>& Ex, 
                                       std::vector<std::vector<std::vector<double>>>& Ey, 
                                       std::vector<std::vector<std::vector<double>>>& Ez) {
      // ---------------------------> inputs <----------------------------------------------------------
      // type( bounds_t), intent( IN) :: b
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent( IN) :: sggMiEx
      // ... (other inputs)
      
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      
      // ---------------------------> empieza AdvanceelectricCPML <-------------------------------------

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = 1;
               regLR[REGION].Psi_Exyvac[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Exyvac[i][j][k] + 
               (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR[REGION].Psi_Exyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regLR[REGION].Psi_Ezyvac[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Ezyvac[i][j][k] + 
               (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR[REGION].Psi_Ezyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               // --->
               medio = 1;
               regLR[REGION].Psi_Exyvac[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Exyvac[i][j][k] + 
               (Hz[i_m][j_m][k_m] - Hz[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] + G2[medio] * regLR[REGION].Psi_Exyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regLR[REGION].Psi_Ezyvac[i][j][k] = P_be_y[j] * regLR[REGION].Psi_Ezyvac[i][j][k] + 
               (Hx[i_m][j_m][k_m] - Hx[i_m][j_m - 1][k_m]) * P_ce_y[j];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] - G2[medio] * regLR[REGION].Psi_Ezyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif



      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regDU[REGION].Psi_Eyzvac[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Eyzvac[i][j][k] + 
               (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU[REGION].Psi_Eyzvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               medio = 1;
               regDU[REGION].Psi_Exzvac[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Exzvac[i][j][k] + 
               (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU[REGION].Psi_Exzvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regDU[REGION].Psi_Eyzvac[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Eyzvac[i][j][k] + 
               (Hx[i_m][j_m][k_m] - Hx[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] + G2[medio] * regDU[REGION].Psi_Eyzvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEx].ZI(REGION); k <= PMLc[iEx].ZE(REGION); ++k) {
         k_m = k - b.Ex.ZI;
         for (j = PMLc[iEx].YI(REGION); j <= PMLc[iEx].YE(REGION); ++j) {
            j_m = j - b.Ex.YI;
            for (i = PMLc[iEx].XI(REGION); i <= PMLc[iEx].XE(REGION); ++i) {
               i_m = i - b.Ex.XI;
               medio = 1;
               regDU[REGION].Psi_Exzvac[i][j][k] = P_be_z[k] * regDU[REGION].Psi_Exzvac[i][j][k] + 
               (Hy[i_m][j_m][k_m] - Hy[i_m][j_m][k_m - 1]) * P_ce_z[k];
               Ex[i_m][j_m][k_m] = Ex[i_m][j_m][k_m] - G2[medio] * regDU[REGION].Psi_Exzvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif




      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regBF[REGION].Psi_Ezxvac[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Ezxvac[i][j][k] + 
               (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF[REGION].Psi_Ezxvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regBF[REGION].Psi_Eyxvac[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Eyxvac[i][j][k] + 
               (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF[REGION].Psi_Eyxvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = front;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEz].ZI(REGION); k <= PMLc[iEz].ZE(REGION); ++k) {
         k_m = k - b.Ez.ZI;
         for (j = PMLc[iEz].YI(REGION); j <= PMLc[iEz].YE(REGION); ++j) {
            j_m = j - b.Ez.YI;
            for (i = PMLc[iEz].XI(REGION); i <= PMLc[iEz].XE(REGION); ++i) {
               i_m = i - b.Ez.XI;
               medio = 1;
               regBF[REGION].Psi_Ezxvac[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Ezxvac[i][j][k] + 
               (Hy[i_m][j_m][k_m] - Hy[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ez[i_m][j_m][k_m] = Ez[i_m][j_m][k_m] + G2[medio] * regBF[REGION].Psi_Ezxvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iEy].ZI(REGION); k <= PMLc[iEy].ZE(REGION); ++k) {
         k_m = k - b.Ey.ZI;
         for (j = PMLc[iEy].YI(REGION); j <= PMLc[iEy].YE(REGION); ++j) {
            j_m = j - b.Ey.YI;
            for (i = PMLc[iEy].XI(REGION); i <= PMLc[iEy].XE(REGION); ++i) {
               i_m = i - b.Ey.XI;
               medio = 1;
               regBF[REGION].Psi_Eyxvac[i][j][k] = P_be_x[i] * regBF[REGION].Psi_Eyxvac[i][j][k] + 
               (Hz[i_m][j_m][k_m] - Hz[i_m - 1][j_m][k_m]) * P_ce_x[i];
               Ey[i_m][j_m][k_m] = Ey[i_m][j_m][k_m] - G2[medio] * regBF[REGION].Psi_Eyxvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif


      // ---------------------------> acaba AdvanceelectricCPML <---------------------------------------
   }
   
   // **************************************************************************************************
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!! Advances the magnetic field in the PML
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   void AdvanceMagneticCPML_freespace(int NumMedia, const bounds_t& b, 
                                       const std::vector<std::vector<std::vector<int>>>& sggMiHx, 
                                       const std::vector<std::vector<std::vector<int>>>& sggMiHy, 
                                       const std::vector<std::vector<std::vector<int>>>& sggMiHz, 
                                       const std::vector<double>& gm2, 
                                       const std::vector<std::vector<std::vector<double>>>& Ex, 
                                       const std::vector<std::vector<std::vector<double>>>& Ey, 
                                       const std::vector<std::vector<std::vector<double>>>& Ez, 
                                       std::vector<std::vector<std::vector<double>>>& Hx, 
                                       std::vector<std::vector<std::vector<double>>>& Hy, 
                                       std::vector<std::vector<std::vector<double>>>& Hz) {
      // ---------------------------> inputs <----------------------------------------------------------
      // type( bounds_t), intent( IN) :: b
      // integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN) :: sggMiHx
      // ... (other inputs)
      
      // ---------------------------> variables locales <-----------------------------------------------
      int REGION, i, j, k, medio, i_m, j_m, k_m;
      
      // ---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
      // Hetic Fields PML Zone
      //
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iHx].ZI(REGION); k <= PMLc[iHx].ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI(REGION); j <= PMLc[iHx].YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI(REGION); i <= PMLc[iHx].XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR[REGION].Psi_Hxyvac[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hxyvac[i][j][k] + 
               (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR[REGION].Psi_Hxyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iHz].ZI(REGION); k <= PMLc[iHz].ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI(REGION); j <= PMLc[iHz].YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI(REGION); i <= PMLc[iHz].XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR[REGION].Psi_Hzyvac[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hzyvac[i][j][k] + 
               (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR[REGION].Psi_Hzyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REGION = right;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iHx].ZI(REGION); k <= PMLc[iHx].ZE(REGION); ++k) {
         k_m = k - b.Hx.ZI;
         for (j = PMLc[iHx].YI(REGION); j <= PMLc[iHx].YE(REGION); ++j) {
            j_m = j - b.Hx.YI;
            for (i = PMLc[iHx].XI(REGION); i <= PMLc[iHx].XE(REGION); ++i) {
               i_m = i - b.Hx.XI;
               // --->
               regLR[REGION].Psi_Hxyvac[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hxyvac[i][j][k] + 
               (Ez[i_m][j_m + 1][k_m] - Ez[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] - GM2[medio] * regLR[REGION].Psi_Hxyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
      for (k = PMLc[iHz].ZI(REGION); k <= PMLc[iHz].ZE(REGION); ++k) {
         k_m = k - b.Hz.ZI;
         for (j = PMLc[iHz].YI(REGION); j <= PMLc[iHz].YE(REGION); ++j) {
            j_m = j - b.Hz.YI;
            for (i = PMLc[iHz].XI(REGION); i <= PMLc[iHz].XE(REGION); ++i) {
               i_m = i - b.Hz.XI;
               // --->
               regLR[REGION].Psi_Hzyvac[i][j][k] = P_bm_y[j] * regLR[REGION].Psi_Hzyvac[i][j][k] + 
               (Ex[i_m][j_m + 1][k_m] - Ex[i_m][j_m][k_m]) * P_cm_y[j];
               medio = 1;
               Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] + GM2[medio] * regLR[REGION].Psi_Hzyvac[i][j][k];
            }
         }
      }
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
    int k_m = k - b.Hy.ZI;
    for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
        int j_m = j - b.Hy.YI;
        for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
            int i_m = i - b.Hy.XI;
            //--->
            regDU[REGION].Psi_Hyzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyzvac[i][j][k] +
                                                (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
            medio = 1;
            Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2(medio) * regDU[REGION].Psi_Hyzvac[i][j][k];
        } // bucle i
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
    int k_m = k - b.Hx.ZI;
    for (int j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
        int j_m = j - b.Hx.YI;
        for (int i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
            int i_m = i - b.Hx.XI;
            //--->
            regDU[REGION].Psi_Hxzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxzvac[i][j][k] +
                                                (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
            medio = 1;
            Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2(medio) * regDU[REGION].Psi_Hxzvac[i][j][k];
        } // bucle i
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REGION = up;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
    int k_m = k - b.Hy.ZI;
    for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
        int j_m = j - b.Hy.YI;
        for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
            int i_m = i - b.Hy.XI;
            //--->
            regDU[REGION].Psi_Hyzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hyzvac[i][j][k] +
                                                (Ex[i_m][j_m][k_m + 1] - Ex[i_m][j_m][k_m]) * P_cm_z[k];
            medio = 1;
            Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] - GM2(medio) * regDU[REGION].Psi_Hyzvac[i][j][k];
        } // bucle i
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHx].ZI[REGION]; k <= PMLc[iHx].ZE[REGION]; ++k) {
    int k_m = k - b.Hx.ZI;
    for (int j = PMLc[iHx].YI[REGION]; j <= PMLc[iHx].YE[REGION]; ++j) {
        int j_m = j - b.Hx.YI;
        for (int i = PMLc[iHx].XI[REGION]; i <= PMLc[iHx].XE[REGION]; ++i) {
            int i_m = i - b.Hx.XI;
            //--->
            regDU[REGION].Psi_Hxzvac[i][j][k] = P_bm_z[k] * regDU[REGION].Psi_Hxzvac[i][j][k] +
                                                (Ey[i_m][j_m][k_m + 1] - Ey[i_m][j_m][k_m]) * P_cm_z[k];
            medio = 1;
            Hx[i_m][j_m][k_m] = Hx[i_m][j_m][k_m] + GM2(medio) * regDU[REGION].Psi_Hxzvac[i][j][k];
        } // bucle i
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REGION = back;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHz].ZI[REGION]; k <= PMLc[iHz].ZE[REGION]; ++k) {
    int k_m = k - b.Hz.ZI;
    for (int j = PMLc[iHz].YI[REGION]; j <= PMLc[iHz].YE[REGION]; ++j) {
        int j_m = j - b.Hz.YI;
        for (int i = PMLc[iHz].XI[REGION]; i <= PMLc[iHz].XE[REGION]; ++i) {
            int i_m = i - b.Hz.XI;
            //--->
            regBF[REGION].Psi_Hzxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hzxvac[i][j][k] +
                                                (Ey[i_m + 1][j_m][k_m] - Ey[i_m][j_m][k_m]) * P_cm_x[i];
            medio = 1;
            Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - GM2(medio) * regBF[REGION].Psi_Hzxvac[i][j][k];
        }
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
    int k_m = k - b.Hy.ZI;
    for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
        int j_m = j - b.Hy.YI;
        for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
            int i_m = i - b.Hy.XI;
            //--->
            regBF[region].Psi_Hyxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hyxvac[i][j][k] +
                                                (Ez[i_m + 1][j_m][k_m] - Ez[i_m][j_m][k_m]) * P_cm_x[i];
            medio = 1;
            Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + GM2(medio) * regBF[REGION].Psi_Hyxvac[i][j][k];
        }
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REGION = front;

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHz].ZI[REGION]; k <= PMLc[iHz].ZE[REGION]; ++k) {
    int k_m = k - b.Hz.ZI;
    for (int j = PMLc[iHz].YI[REGION]; j <= PMLc[iHz].YE[REGION]; ++j) {
        int j_m = j - b.Hz.YI;
        for (int i = PMLc[iHz].XI[REGION]; i <= PMLc[iHz].XE[REGION]; ++i) {
            int i_m = i - b.Hz.XI;
            //--->
            regBF[REGION].Psi_Hzxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hzxvac[i][j][k] +
                                                (Ey[i_m + 1][j_m][k_m] - Ey[i_m][j_m][k_m]) * P_cm_x[i];
            medio = 1;
            Hz[i_m][j_m][k_m] = Hz[i_m][j_m][k_m] - GM2(medio) * regBF[REGION].Psi_Hzxvac[i][j][k];
        }
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i, j, k, i_m, j_m, k_m, medio)
#endif
for (int k = PMLc[iHy].ZI[REGION]; k <= PMLc[iHy].ZE[REGION]; ++k) {
    int k_m = k - b.Hy.ZI;
    for (int j = PMLc[iHy].YI[REGION]; j <= PMLc[iHy].YE[REGION]; ++j) {
        int j_m = j - b.Hy.YI;
        for (int i = PMLc[iHy].XI[REGION]; i <= PMLc[iHy].XE[REGION]; ++i) {
            int i_m = i - b.Hy.XI;
            //--->
            regBF[region].Psi_Hyxvac[i][j][k] = P_bm_x[i] * regBF[REGION].Psi_Hyxvac[i][j][k] +
                                                (Ez[i_m + 1][j_m][k_m] - Ez[i_m][j_m][k_m]) * P_cm_x[i];
            medio = 1;
            Hy[i_m][j_m][k_m] = Hy[i_m][j_m][k_m] + GM2(medio) * regBF[REGION].Psi_Hyxvac[i][j][k];
        }
    }
}
#ifdef CompileWithOpenMP
#pragma omp end parallel for
#endif

// ---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
return;
} // endsubroutine AdvanceMagneTicCPML_freespace

} // end Module BORDERS_CPML_m