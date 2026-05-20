#include <vector>
#include <complex>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

using RKIND = double;
using CKIND = std::complex<double>;
using INTEGERSIZEOFMEDIAMATRICES = int;
using RKIND_tiempo = double;
constexpr int BUFSIZE = 256;
constexpr int iEx = 1, iEy = 2, iEz = 3, iHx = 4, iHy = 5, iHz = 6;
constexpr double mcpi2 = 6.28318530717958647693;
constexpr double pi = 3.14159265358979323846;
constexpr int REALSIZE = 8;

namespace External {
    struct bounds_t {
        struct { int XI, XE, YI, YE, ZI, ZE, NX, NY, NZ; } Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh;
    };
    struct SGGFDTDINFO_t {
        struct { int XI, XE, YI, YE, ZI, ZE; } alloc[7];
        RKIND dt;
        std::vector<RKIND> LineX, LineY, LineZ;
        struct { int XI, XE, YI, YE, ZI, ZE; } SINPMLSweep[7];
        struct { bool IsBackPEC, IsBackPMC, IsFrontPEC, IsFrontPMC, IsLeftPEC, IsLeftPMC, IsRightPEC, IsRightPMC, IsDownPEC, IsDownPMC, IsUpPEC, IsUpPMC; } Border;
        struct { bool is_pec; } med[100];
    };
    struct limit_t { int XI, XE, YI, YE, ZI, ZE; };
    struct nf2ff_t { bool Tr, Fr, Iz, De, Ab, Ar; };
    struct coorsxyzP_t {
        struct { std::vector<RKIND>* x; std::vector<RKIND>* y; std::vector<RKIND>* z; } PhysCoor[7];
    };
    void stoponerror(int, int, const std::string&) {}
    void print11(int, const std::string&, bool = false) {}
    const std::string SEPARADOR = "================================";
}

namespace farfield_m {

    struct ehxyz_t {
        int Ex = -15, Ey = -15, Ez = -15, Hx = -15, Hy = -15, Hz = -15;
    };

    struct tfidaa_t {
        ehxyz_t com, fin, tra, fro, izq, der, aba, arr;
    };

    struct ijk_t {
        tfidaa_t i, j, k;
    };

    struct co_t {
        RKIND x_Mx, y_Mx, z_Mx;
        RKIND x_My, y_My, z_My;
        RKIND x_Mz, y_Mz, z_Mz;
        RKIND x_Jx, y_Jx, z_Jx;
        RKIND x_Jy, y_Jy, z_Jy;
        RKIND x_Jz, y_Jz, z_Jz;
    };

    struct farfield_t {
        ijk_t TrFr, IzDe, AbAr;
        bool farfieldTr = false, farfieldFr = false, farfieldIz = false, farfieldDe = false, farfieldAr = false, farfieldAb = false;
        bool farfieldTr_ClonePEC_Front = false, farfieldTr_ClonePEC_Left = false, farfieldTr_ClonePEC_Right = false, farfieldTr_ClonePEC_Up = false, farfieldTr_ClonePEC_Down = false;
        bool farfieldTr_ClonePMC_Front = false, farfieldTr_ClonePMC_Left = false, farfieldTr_ClonePMC_Right = false, farfieldTr_ClonePMC_Up = false, farfieldTr_ClonePMC_Down = false;
        bool farfieldFr_ClonePEC_Back = false, farfieldFr_ClonePEC_Left = false, farfieldFr_ClonePEC_Right = false, farfieldFr_ClonePEC_Up = false, farfieldFr_ClonePEC_Down = false;
        bool farfieldFr_ClonePMC_Back = false, farfieldFr_ClonePMC_Left = false, farfieldFr_ClonePMC_Right = false, farfieldFr_ClonePMC_Up = false, farfieldFr_ClonePMC_Down = false;
        bool farfieldIz_ClonePEC_Back = false, farfieldIz_ClonePEC_Front = false, farfieldIz_ClonePEC_Right = false, farfieldIz_ClonePEC_Up = false, farfieldIz_ClonePEC_Down = false;
        bool farfieldIz_ClonePMC_Back = false, farfieldIz_ClonePMC_Front = false, farfieldIz_ClonePMC_Right = false, farfieldIz_ClonePMC_Up = false, farfieldIz_ClonePMC_Down = false;
        bool farfieldDe_ClonePEC_Back = false, farfieldDe_ClonePEC_Front = false, farfieldDe_ClonePEC_Left = false, farfieldDe_ClonePEC_Up = false, farfieldDe_ClonePEC_Down = false;
        bool farfieldDe_ClonePMC_Back = false, farfieldDe_ClonePMC_Front = false, farfieldDe_ClonePMC_Left = false, farfieldDe_ClonePMC_Up = false, farfieldDe_ClonePMC_Down = false;
        bool farfieldAr_ClonePEC_Back = false, farfieldAr_ClonePEC_Front = false, farfieldAr_ClonePEC_Left = false, farfieldAr_ClonePEC_Right = false, farfieldAr_ClonePEC_Down = false;
        bool farfieldAr_ClonePMC_Back = false, farfieldAr_ClonePMC_Front = false, farfieldAr_ClonePMC_Left = false, farfieldAr_ClonePMC_Right = false, farfieldAr_ClonePMC_Down = false;
        bool farfieldAb_ClonePEC_Back = false, farfieldAb_ClonePEC_Front = false, farfieldAb_ClonePEC_Left = false, farfieldAb_ClonePEC_Right = false, farfieldAb_ClonePEC_Up = false;
        bool farfieldAb_ClonePMC_Back = false, farfieldAb_ClonePMC_Front = false, farfieldAb_ClonePMC_Left = false, farfieldAb_ClonePMC_Right = false, farfieldAb_ClonePMC_Up = false;
        std::vector<std::vector<std::vector<CKIND>>> ExIz, ExDe, ExAb, ExAr, EyFr, EyTr, EyAb, EyAr, EzIz, EzDe, EzFr, EzTr;
        std::vector<std::vector<std::vector<CKIND>>> HxIz, HxDe, HxAb, HxAr, HyFr, HyTr, HyAb, HyAr, HzIz, HzDe, HzFr, HzTr;
        std::vector<std::vector<std::vector<CKIND>>> HxIz2, HxDe2, HxAb2, HxAr2, HyFr2, HyTr2, HyAb2, HyAr2, HzIz2, HzDe2, HzFr2, HzTr2;
        std::vector<CKIND> expIwdt, auxExp_E, auxExp_H, dftEntrada;
        int NumFreqs = 0, esqx1 = 0, esqx2 = 0, esqy1 = 0, esqy2 = 0, esqz1 = 0, esqz2 = 0, Ndecim = 0;
        External::coorsxyzP_t Punto;
        RKIND InitialFreq = 0, FinalFreq = 0, FreqStep = 0, dtDecim = 0;
        RKIND thetaStart = 0, thetaStop = 0, thetaStep = 0;
        RKIND phiStart = 0, phiStop = 0, phiStep = 0;
        std::string FileNormalize;
        int unitfarfield = 0;
        std::string filefarfield;
        RKIND XDobleAncho = 0, YDobleAncho = 0, ZDobleAncho = 0;
        RKIND XOffsetPlus = 0, YOffsetPlus = 0, ZOffsetPlus = 0;
        RKIND XOffsetMinus = 0, YOffsetMinus = 0, ZOffsetMinus = 0;
#ifdef CompileWithMPI
        int MPISubComm = 0, MPIRoot = 0;
#endif
    };

    RKIND cluz = 0, zvac = 0;
    RKIND eps0 = 0, mu0 = 0;
    farfield_t FF;

    CKIND average(int, const CKIND&, const CKIND&) { return CKIND(0.0, 0.0); }

    void InitFarField(
        const External::SGGFDTDINFO_t&,
        const std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>>&,
        const std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>>&,
        const std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>>&,
        const std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>>&,
        const std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>>&,
        const std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>>&,
        int, int, const External::bounds_t&, bool, int, const std::string&,
        int, int, int, int, int, int,
        RKIND, RKIND, RKIND,
        RKIND, RKIND, RKIND,
        RKIND, RKIND, RKIND,
        const std::string&,
        const std::vector<External::limit_t>&,
        const External::nf2ff_t&, bool, RKIND, RKIND
#ifdef CompileWithMPI
        , int, int
#endif
    ) {}

    void UpdateFarField(int, const External::bounds_t&,
        const std::vector<std::vector<std::vector<RKIND>>>&,
        const std::vector<std::vector<std::vector<RKIND>>>&,
        const std::vector<std::vector<std::vector<RKIND>>>&,
        const std::vector<std::vector<std::vector<RKIND>>>&,
        const std::vector<std::vector<std::vector<RKIND>>>&,
        const std::vector<std::vector<std::vector<RKIND>>>&) {}

}
if (j <= FF.ABAr.J.fin.Ey)  Mx = - EcampoY(i_m, j_m, ii) * dxh(i_m) * dye(j_m) * NORMAL;
                              if (i <= FF.ABAr.I.fin.Ex)  My = + EcampoX(i_m, j_m, ii) * dxe(i_m) * dyh(j_m) * NORMAL;
                              if (j <= FF.ABAr.J.fin.Hy)  Jx = + Average(pasadas, HcampoY(i_m, j_m, ii), Hcampo2Y(i_m, j_m, ii)) * dxe(i_m) * dyh(j_m) * NORMAL;
                              if (i <= FF.ABAr.I.fin.Hx)  Jy = - Average(pasadas, HcampoX(i_m, j_m, ii), Hcampo2X(i_m, j_m, ii)) * dxh(i_m) * dye(j_m) * NORMAL;
                              co.x_Mx = FF.Punto.PhysCoor(iEy).x(i); co.y_Mx = FF.Punto.PhysCoor(iEy).y(j); co.z_Mx = FF.Punto.PhysCoor(iEy).z(k);
                              co.x_My = FF.Punto.PhysCoor(iEx).x(i); co.y_My = FF.Punto.PhysCoor(iEx).y(j); co.z_My = FF.Punto.PhysCoor(iEx).z(k);
                              co.x_Jx = co.x_My;                     co.y_Jx = co.y_My;                     co.z_Jx = co.z_My;
                              co.x_Jy = co.x_Mx;                     co.y_Jy = co.y_Mx;                     co.z_Jy = co.z_Mx;
                              update_LN(comun, co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, Mx, My, Mz, Jx, Jy, Jz, L_theta, L_phi, N_theta, N_phi);
                              // simetrias
                              // la AbAr hay que llamarla para cada caso
                              new_Mx = Mx;
                              new_My = My;
                              new_Jx = Jx;
                              new_Jy = Jy;
                              new_co = co;
                              cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              //
                              if (FF.farfieldAb_ClonePEC_LEFT || FF.farfieldAr_ClonePEC_LEFT) {
                                 // ?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
                                 new_Mx = + Mx;
                                 new_My = - My;
                                 new_Jx = - Jx;
                                 new_Jy = + Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetMinus;
                                 new_co.y_My = -co.y_My + FF.yOffsetMinus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              if (FF.farfieldAb_ClonePMC_LEFT || FF.farfieldAr_ClonePMC_LEFT) {
                                 new_Mx = - Mx;
                                 new_My = + My;
                                 new_Jx = + Jx;
                                 new_Jy = - Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetMinus;
                                 new_co.y_My = -co.y_My + FF.yOffsetMinus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              //
                              if (FF.farfieldAb_ClonePEC_RIGHT || FF.farfieldAr_ClonePEC_RIGHT) {
                                 new_Mx = + Mx;
                                 new_My = - My;
                                 new_Jx = - Jx;
                                 new_Jy = + Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetPlus;
                                 new_co.y_My = -co.y_My + FF.yOffsetPlus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              if (FF.farfieldAb_ClonePMC_RIGHT || FF.farfieldAr_ClonePMC_RIGHT) {
                                 new_Mx = - Mx;
                                 new_My = + My;
                                 new_Jx = + Jx;
                                 new_Jy = - Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetPlus;
                                 new_co.y_My = -co.y_My + FF.yOffsetPlus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              //
                              if (FF.farfieldAb_ClonePEC_BACK || FF.farfieldAr_ClonePEC_BACK) {
                                 new_Mx = - Mx;
                                 new_My = + My;
                                 new_Jx = + Jx;
                                 new_Jy = - Jy;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetMinus;
                                 new_co.x_My = -co.x_My + FF.xOffsetMinus;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              if (FF.farfieldAb_ClonePMC_BACK || FF.farfieldAr_ClonePMC_BACK) {
                                 new_Mx = + Mx;
                                 new_My = - My;
                                 new_Jx = - Jx;
                                 new_Jy = + Jy;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetMinus;
                                 new_co.x_My = -co.x_My + FF.xOffsetMinus;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              //
                              if (FF.farfieldAb_ClonePEC_FRONT || FF.farfieldAr_ClonePEC_FRONT) {
                                 new_Mx = - Mx;
                                 new_My = + My;
                                 new_Jx = + Jx;
                                 new_Jy = - Jy;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetPlus;
                                 new_co.x_My = -co.x_My + FF.xOffsetPlus;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              if (FF.farfieldAb_ClonePMC_FRONT || FF.farfieldAr_ClonePMC_FRONT) {
                                 new_Mx = + Mx;
                                 new_My = - My;
                                 new_Jx = - Jx;
                                 new_Jy = + Jy;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetPlus;
                                 new_co.x_My = -co.x_My + FF.xOffsetPlus;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }

                              // CASOS MIXTOS esquinas
                              if (((FF.farfieldAb_ClonePEC_LEFT || FF.farfieldAr_ClonePEC_LEFT) && (FF.farfieldAb_ClonePEC_BACK || FF.farfieldAr_ClonePEC_BACK)) ||
                                  ((FF.farfieldAb_ClonePMC_LEFT || FF.farfieldAr_ClonePMC_LEFT) && (FF.farfieldAb_ClonePMC_BACK || FF.farfieldAr_ClonePMC_BACK)) ||
                                  ((FF.farfieldAb_ClonePEC_LEFT || FF.farfieldAr_ClonePEC_LEFT) && (FF.farfieldAb_ClonePMC_BACK || FF.farfieldAr_ClonePMC_BACK)) ||
                                  ((FF.farfieldAb_ClonePMC_LEFT || FF.farfieldAr_ClonePMC_LEFT) && (FF.farfieldAb_ClonePEC_BACK || FF.farfieldAr_ClonePEC_BACK))) {
                                 sigNo = +1.0_RKIND;
                                 if (((FF.farfieldAb_ClonePEC_LEFT || FF.farfieldAr_ClonePEC_LEFT) && (FF.farfieldAb_ClonePEC_BACK || FF.farfieldAr_ClonePEC_BACK)) ||
                                     ((FF.farfieldAb_ClonePMC_LEFT || FF.farfieldAr_ClonePMC_LEFT) && (FF.farfieldAb_ClonePMC_BACK || FF.farfieldAr_ClonePMC_BACK))) {
                                    sigNo = -1.0_RKIND;
                                 }
                                 new_Mx = signo * Mx;
                                 new_My = signo * My;
                                 new_Jx = signo * Jx;
                                 new_Jy = signo * Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetMinus;
                                 new_co.y_My = -co.y_My + FF.yOffsetMinus;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetMinus;
                                 new_co.x_My = -co.x_My + FF.xOffsetMinus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }

                              if (((FF.farfieldAb_ClonePEC_LEFT || FF.farfieldAr_ClonePEC_LEFT) && (FF.farfieldAb_ClonePEC_FRONT || FF.farfieldAr_ClonePEC_FRONT)) ||
                                  ((FF.farfieldAb_ClonePMC_LEFT || FF.farfieldAr_ClonePMC_LEFT) && (FF.farfieldAb_ClonePMC_FRONT || FF.farfieldAr_ClonePMC_FRONT)) ||
                                  ((FF.farfieldAb_ClonePEC_LEFT || FF.farfieldAr_ClonePEC_LEFT) && (FF.farfieldAb_ClonePMC_FRONT || FF.farfieldAr_ClonePMC_FRONT)) ||
                                  ((FF.farfieldAb_ClonePMC_LEFT || FF.farfieldAr_ClonePMC_LEFT) && (FF.farfieldAb_ClonePEC_FRONT || FF.farfieldAr_ClonePEC_FRONT))) {
                                 sigNo = +1.0_RKIND;
                                 if (((FF.farfieldAb_ClonePEC_LEFT || FF.farfieldAr_ClonePEC_LEFT) && (FF.farfieldAb_ClonePEC_FRONT || FF.farfieldAr_ClonePEC_FRONT)) ||
                                     ((FF.farfieldAb_ClonePMC_LEFT || FF.farfieldAr_ClonePMC_LEFT) && (FF.farfieldAb_ClonePMC_FRONT || FF.farfieldAr_ClonePMC_FRONT))) {
                                    sigNo = -1.0_RKIND;
                                 }
                                 new_Mx = signo * Mx;
                                 new_My = signo * My;
                                 new_Jx = signo * Jx;
                                 new_Jy = signo * Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetMinus;
                                 new_co.y_My = -co.y_My + FF.yOffsetMinus;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetPlus;
                                 new_co.x_My = -co.x_My + FF.xOffsetPlus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);

cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }

                              if (((FF.farfieldAb_ClonePEC_RIGHT || FF.farfieldAr_ClonePEC_RIGHT) && (FF.farfieldAb_ClonePEC_BACK || FF.farfieldAr_ClonePEC_BACK)) ||
                                  ((FF.farfieldAb_ClonePMC_RIGHT || FF.farfieldAr_ClonePMC_RIGHT) && (FF.farfieldAb_ClonePMC_BACK || FF.farfieldAr_ClonePMC_BACK)) ||
                                  ((FF.farfieldAb_ClonePEC_RIGHT || FF.farfieldAr_ClonePEC_RIGHT) && (FF.farfieldAb_ClonePMC_BACK || FF.farfieldAr_ClonePMC_BACK)) ||
                                  ((FF.farfieldAb_ClonePMC_RIGHT || FF.farfieldAr_ClonePMC_RIGHT) && (FF.farfieldAb_ClonePEC_BACK || FF.farfieldAr_ClonePEC_BACK))) {
                                 sigNo = +1.0_RKIND;
                                 if (((FF.farfieldAb_ClonePEC_RIGHT || FF.farfieldAr_ClonePEC_RIGHT) && (FF.farfieldAb_ClonePEC_BACK || FF.farfieldAr_ClonePEC_BACK)) ||
                                     ((FF.farfieldAb_ClonePMC_RIGHT || FF.farfieldAr_ClonePMC_RIGHT) && (FF.farfieldAb_ClonePMC_BACK || FF.farfieldAr_ClonePMC_BACK))) {
                                    sigNo = -1.0_RKIND;
                                 }
                                 new_Mx = signo * Mx;
                                 new_My = signo * My;
                                 new_Jx = signo * Jx;
                                 new_Jy = signo * Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetPlus;
                                 new_co.y_My = -co.y_My + FF.yOffsetPlus;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetMinus;
                                 new_co.x_My = -co.x_My + FF.xOffsetMinus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }

                              if (((FF.farfieldAb_ClonePEC_RIGHT || FF.farfieldAr_ClonePEC_RIGHT) && (FF.farfieldAb_ClonePEC_FRONT || FF.farfieldAr_ClonePEC_FRONT)) ||
                                  ((FF.farfieldAb_ClonePMC_RIGHT || FF.farfieldAr_ClonePMC_RIGHT) && (FF.farfieldAb_ClonePMC_FRONT || FF.farfieldAr_ClonePMC_FRONT)) ||
                                  ((FF.farfieldAb_ClonePEC_RIGHT || FF.farfieldAr_ClonePEC_RIGHT) && (FF.farfieldAb_ClonePMC_FRONT || FF.farfieldAr_ClonePMC_FRONT)) ||
                                  ((FF.farfieldAb_ClonePMC_RIGHT || FF.farfieldAr_ClonePMC_RIGHT) && (FF.farfieldAb_ClonePEC_FRONT || FF.farfieldAr_ClonePEC_FRONT))) {
                                 sigNo = +1.0_RKIND;
                                 if (((FF.farfieldAb_ClonePEC_RIGHT || FF.farfieldAr_ClonePEC_RIGHT) && (FF.farfieldAb_ClonePEC_FRONT || FF.farfieldAr_ClonePEC_FRONT)) ||
                                     ((FF.farfieldAb_ClonePMC_RIGHT || FF.farfieldAr_ClonePMC_RIGHT) && (FF.farfieldAb_ClonePMC_FRONT || FF.farfieldAr_ClonePMC_FRONT))) {
                                    sigNo = -1.0_RKIND;
                                 }
                                 new_Mx = signo * Mx;
                                 new_My = signo * My;
                                 new_Jx = signo * Jx;
                                 new_Jy = signo * Jy;
                                 new_co.y_Mx = -co.y_Mx + FF.yOffsetPlus;
                                 new_co.y_My = -co.y_My + FF.yOffsetPlus;
                                 new_co.x_Mx = -co.x_Mx + FF.xOffsetPlus;
                                 new_co.x_My = -co.x_My + FF.xOffsetPlus;
                                 new_co.y_Jx = new_co.y_My;
                                 new_co.y_Jy = new_co.y_Mx;
                                 new_co.x_Jx = new_co.x_My;
                                 new_co.x_Jy = new_co.x_Mx;
                                 update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
                                 cloneAbAr(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi, NORMAL);
                              }
                              // end symmetries
                           }
                        }
                        L_theta_final = L_theta_final + L_theta; L_phi_final = L_phi_final + L_phi;
                        N_theta_final = N_theta_final + N_theta; N_phi_final = N_phi_final + N_phi;
                        L_theta = 0.0_RKIND; L_phi = 0.0_RKIND; N_theta = 0.0_RKIND; N_phi = 0.0_RKIND;
                     } // end goAhead
                  } // end donde

                  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CompileWithMPI
                  MPI_Barrier(FF.MPISubComm, ierr);
                  dummy = real(L_theta_final);
                  MPI_AllReduce(&dummy, &newdummy1, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);
                  dummy = AIMAG(L_theta_final);
                  MPI_AllReduce(&dummy, &newdummy2, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);

                  L_theta_final = newdummy1 + (0.0_RKIND, 1.0_RKIND) * newdummy2;
                  //
                  dummy = real(L_phi_final);
                  MPI_AllReduce(&dummy, &newdummy1, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);
                  dummy = AIMAG(L_phi_final);
                  MPI_AllReduce(&dummy, &newdummy2, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);
                  //
                  L_phi_final = newdummy1 + (0.0_RKIND, 1.0_RKIND) * newdummy2;

                  dummy = real(N_theta_final);
                  MPI_AllReduce(&dummy, &newdummy1, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);
                  dummy = AIMAG(N_theta_final);
                  MPI_AllReduce(&dummy, &newdummy2, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);

                  N_theta_final = newdummy1 + (0.0_RKIND, 1.0_RKIND) * newdummy2;
                  //
                  dummy = real(N_phi_final);
                  MPI_AllReduce(&dummy, &newdummy1, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);
                  dummy = AIMAG(N_phi_final);
                  MPI_AllReduce(&dummy, &newdummy2, 1, REALSIZE, MPI_SUM, FF.MPISubComm, &ierr);
                  MPI_Barrier(FF.MPISubComm, &ierr);
                  //
                  N_phi_final = newdummy1 + (0.0_RKIND, 1.0_RKIND) * newdummy2;
                  MPI_Barrier(FF.MPISubComm, &ierr);
#endif
                  Etheta(pasadas) = -(0, 1.0_RKIND) * freq / (2.0_RKIND * cluz) * (L_phi_final + zvac * N_theta_final); // /FF.dftEntrada(ii) !no normalize to calculate power correctly
                  Ephi(pasadas) = (0, 1.0_RKIND) * freq / (2.0_RKIND * cluz) * (L_theta_final - zvac * N_phi_final); // /FF.dftEntrada(ii) !no normalize to calculate power correctly
                  RCS(pasadas) = (2.0_RKIND * pi * freq / cluz) ** 2.0_RKIND / (4.0_RKIND * pi * Abs(FF.dftEntrada(ii)) ** 2.0_RKIND) *
                                 (Abs(L_phi_final + zvac * N_theta_final) ** 2.0_RKIND + Abs(L_theta_final - zvac * N_phi_final) ** 2.0_RKIND);


#ifdef CompileWithMPI
                  if (FF.MPIRoot == layoutnumber) {
#endif
                  if (pasadas == 1) {
                     fprintf(FF.unitfarfield, fmt, freq, theta, phi,
                             Abs(Etheta(2)), ATAN2(AIMAG(Etheta(2)), real(Etheta(2))), // PASADAS=2=GEOMETRIC, PASADAS=1=ARITHMETIC
                             Abs(Ephi(2)), ATAN2(AIMAG(Ephi(2)), real(Ephi(2))), RCS(1), RCS(2));
                  }

#ifdef CompileWithMPI
                  }
#endif
               } // end pasadas
            } // end phi loop
         } // end theta loop

         sprintf(dubuf, "NF2FF: End processing freq= %e", freq);
         // call print11(layoutnumber, dubuf, .TRUE.)

      } // end frequency loop
      // end calculation
      if (Thetavector.is_allocated()) Thetavector.deallocate();
      if (sizephi.is_allocated()) sizephi.deallocate();
      if (Phimatrix.is_allocated()) Phimatrix.deallocate();

      fclose(FF.unitfarfield);



      sprintf(dubuf, "NF2FF: END ");
      if (layoutnumber == 0) print11(layoutnumber, dubuf, true);

   } // end subroutine


   void update_LN(co_t co, double sintheta_cosphi, double sintheta_sinphi, double costheta, double costheta_cosphi, double costheta_sinphi, double sintheta, double sinphi, double cosphi, complex<double> Mx, complex<double> My, complex<double> Mz, complex<double> Jx, complex<double> Jy, complex<double> Jz, complex<double>& L_theta, complex<double>& L_phi, complex<double>& N_theta, complex<double>& N_phi) {

      complex<double> comunMx, comunMy, comunMz, comunJx, comunJy, comunJz, comun;
      // might be missing a sign in the exponential although it doesn't affect anything (I think it's -exp[-j k r] taflove 3rd 372 nf2ff) 02/03/15
      comunMx = exp(comun * (co.x_Mx * sintheta_cosphi + co.y_Mx * sintheta_sinphi + co.z_Mx * costheta));
      comunMy = exp(comun * (co.x_My * sintheta_cosphi + co.y_My * sintheta_sinphi + co.z_My * costheta));
      comunMz = exp(comun * (co.x_Mz * sintheta_cosphi + co.y_Mz * sintheta_sinphi + co.z_Mz * costheta));
      comunJx = exp(comun * (co.x_Jx * sintheta_cosphi + co.y_Jx * sintheta_sinphi + co.z_Jx * costheta));
      comunJy = exp(comun * (co.x_Jy * sintheta_cosphi + co.y_Jy * sintheta_sinphi + co.z_Jy * costheta));
      comunJz = exp(comun * (co.x_Jz * sintheta_cosphi + co.y_Jz * sintheta_sinphi + co.z_Jz * costheta));
      //
      L_theta = L_theta + (Mx * costheta_cosphi * comunMx + My * costheta_sinphi * comunMy - Mz * sintheta * comunMz);
      L_phi = L_phi + (-Mx * sinphi * comunMx + My * cosphi * comunMy);
      //
      N_theta = N_theta + (Jx * costheta_cosphi * comunJx + Jy * costheta_sinphi * comunJy - Jz * sintheta * comunJz);
      N_phi = N_phi + (-Jx * sinphi * comunJx + Jy * cosphi * comunJy);
   }



   void cloneTrFr(co_t co, double sintheta_cosphi, double sintheta_sinphi, double costheta, double costheta_cosphi, double costheta_sinphi, double sintheta, double sinphi, double cosphi, complex<double> Mx, complex<double> My, complex<double> Mz, complex<double> Jx, complex<double> Jy, complex<double> Jz, complex<double>& L_theta, complex<double>& L_phi, complex<double>& N_theta, complex<double>& N_phi, double NORMAL) {
      co_t new_co;
      complex<double> new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz;
      complex<double> comun;

      new_co = co; new_Mx = Mx; new_My = My; new_Mz = Mz; new_Jx = Jx; new_Jy = Jy; new_Jz = Jz;
      if (FF.farfieldTr_ClonePEC_Front || FF.farfieldFr_ClonePEC_Back) {
         new_My = +My; // only in this case normals change
         new_Mz = +Mz;
         new_Jy = -Jy;
         new_Jz = -Jz;
         new_co.x_My = co.x_My + FF.XDobleAncho * NORMAL; // sign change add or subtract distance
         new_co.x_Mz = co.x_Mz + FF.XDobleAncho * NORMAL;
         new_co.x_Jy = new_co.x_Mz;
         new_co.x_Jz = new_co.x_My;
         update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
      }
      if (FF.farfieldTr_ClonePMC_Front || FF.farfieldFr_ClonePMC_Back) {
         new_My = -My;
         new_Mz = -Mz;
         new_Jy = +Jy;
         new_Jz = +Jz;
         new_co.x_My = co.x_My + FF.XDobleAncho * NORMAL; // sign change add or subtract distance
         new_co.x_Mz = co.x_Mz + FF.XDobleAncho * NORMAL;
         new_co.x_Jy = new_co.x_Mz;
         new_co.x_Jz = new_co.x_My;
         update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
      }
      return;
   }



   void cloneIzDe(co_t co, double sintheta_cosphi, double sintheta_sinphi, double costheta, double costheta_cosphi, double costheta_sinphi, double sintheta, double sinphi, double cosphi, complex<double> Mx, complex<double> My, complex<double> Mz, complex<double> Jx, complex<double> Jy, complex<double> Jz, complex<double>& L_theta, complex<double>& L_phi, complex<double>& N_theta, complex<double>& N_phi, double NORMAL) {
      co_t new_co;
      complex<double> new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz;
      complex<double> comun;

      new_co = co; new_Mx = Mx; new_My = My; new_Mz = Mz; new_Jx = Jx; new_Jy = Jy; new_Jz = Jz;
      if (FF.farfieldIz_ClonePEC_Right || FF.farfieldDe_ClonePEC_Left) {
         new_Mx = +Mx; // only in this case normals change
         new_Mz = +Mz;
         new_Jx = -Jx;
         new_Jz = -Jz;
         new_co.y_Mx = co.y_Mx + FF.YDobleAncho * NORMAL; // sign change add or subtract distance
         new_co.y_Mz = co.y_Mz + FF.YDobleAncho * NORMAL;
         new_co.y_Jx = new_co.y_Mz;
         new_co.y_Jz = new_co.y_Mx;
         update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
      }
      if (FF.farfieldIz_ClonePMC_Right || FF.farfieldDe_ClonePMC_Left) {
         new_Mx = -Mx; // only in this case normals change
         new_Mz = -Mz;
         new_Jx = +Jx;
         new_Jz = +Jz;
         new_co.y_Mx = co.y_Mx + FF.YDobleAncho * NORMAL; // sign change add or subtract distance
         new_co.y_Mz = co.y_Mz + FF.YDobleAncho * NORMAL;
         new_co.y_Jx = new_co.y_Mz;
         new_co.y_Jz = new_co.y_Mx;
         update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
      }
      return;
   }


   void cloneAbAr(co_t co, double sintheta_cosphi, double sintheta_sinphi, double costheta, double costheta_cosphi, double costheta_sinphi, double sintheta, double sinphi, double cosphi, complex<double> Mx, complex<double> My, complex<double> Mz, complex<double> Jx, complex<double> Jy, complex<double> Jz, complex<double>& L_theta, complex<double>& L_phi, complex<double>& N_theta, complex<double>& N_phi, double NORMAL) {
      co_t new_co;
      complex<double> new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz;
      complex<double> comun;

#include <complex>
#include <cmath>
#include <vector>

// Assuming necessary types and constants are defined in headers or elsewhere
// RKIND, CKIND, PI, FF, comun, update_LN, etc. need to be resolved from context.
// For the purpose of this translation, we assume:
// - RKIND is double
// - CKIND is std::complex<double>
// - FF is an object/struct with the specified members
// - comun is passed by reference or is a global/context object
// - update_LN is a function
// - co, Mx, My, Mz, Jx, Jy, Jz, L_theta, L_phi, N_theta, N_phi are variables/objects

// Forward declarations or assumptions based on typical usage in such modules
// struct CommonData; // comun
// struct CoordType;  // co, new_co
// struct FarFieldData; // FF
// void update_LN(CommonData& comun, CoordType& new_co, double& sintheta_cosphi, double& sintheta_sinphi, double& costheta, double& costheta_cosphi, double& costheta_sinphi, double& sintheta, double& sinphi, double& cosphi, double& new_Mx, double& new_My, double& new_Mz, double& new_Jx, double& new_Jy, double& new_Jz, std::vector<double>& L_theta, std::vector<double>& L_phi, std::vector<double>& N_theta, std::vector<double>& N_phi);

namespace farfield_m {

void some_subroutine(
    // Assuming these are inputs/outputs to the subroutine containing the first block
    // The original code snippet starts mid-subroutine, so we define a placeholder
    // Inputs from context not fully visible, but used in assignments
    const double& Mx, const double& My, const double& Mz,
    const double& Jx, const double& Jy, const double& Jz,
    CoordType& co, // Assuming co is passed by reference or is a member
    FarFieldData& FF,
    CommonData& comun,
    std::vector<double>& L_theta,
    std::vector<double>& L_phi,
    std::vector<double>& N_theta,
    std::vector<double>& N_phi,
    double& NORMAL
) {
    double sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi;
    double new_Mx = Mx;
    double new_My = My;
    double new_Mz = Mz;
    double new_Jx = Jx;
    double new_Jy = Jy;
    double new_Jz = Jz;
    CoordType new_co = co; // Assuming copy assignment is valid for CoordType

    if (FF.farfieldAb_ClonePEC_UP || FF.farfieldAr_ClonePEC_DOWN) {
        new_Mx = +Mx;
        new_My = +My;
        new_Jx = -Jx;
        new_Jy = -Jy;
        new_co.Z_Mx = co.Z_Mx + FF.ZDobleAncho * NORMAL;
        new_co.Z_My = co.Z_My + FF.ZDobleAncho * NORMAL;
        new_co.Z_Jx = new_co.Z_My;
        new_co.Z_Jy = new_co.Z_Mx;
        update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
    }

    if (FF.farfieldAb_ClonePMC_UP || FF.farfieldAr_ClonePMC_DOWN) {
        new_Mx = -Mx;
        new_My = -My;
        new_Jx = +Jx;
        new_Jy = +Jy;
        new_co.Z_Mx = co.Z_Mx + FF.ZDobleAncho * NORMAL;
        new_co.Z_My = co.Z_My + FF.ZDobleAncho * NORMAL;
        new_co.Z_Jx = new_co.Z_My;
        new_co.Z_Jy = new_co.Z_Mx;
        update_LN(comun, new_co, sintheta_cosphi, sintheta_sinphi, costheta, costheta_cosphi, costheta_sinphi, sintheta, sinphi, cosphi, new_Mx, new_My, new_Mz, new_Jx, new_Jy, new_Jz, L_theta, L_phi, N_theta, N_phi);
    }
}

std::complex<double> average(int pasadas, std::complex<double> z1, std::complex<double> z2) {
    std::complex<double> Z(0.0, 0.0);
    double phi1, phi2, nphi1, nphi2;

    if (pasadas == 2) { // geometrica
        phi1 = std::atan2(std::imag(z1), std::real(z1));
        phi2 = std::atan2(std::imag(z2), std::real(z2));

        // TRAMPA CHINA PARA EVITAR LOS BRANCH CUT EN 0 Y PI
        nphi1 = phi1;
        nphi2 = phi2;
        if ((phi1 < -M_PI / 2.0) && (phi2 > M_PI / 2.0)) {
            nphi1 = phi1 - 2.0 * M_PI;
        }
        if ((phi2 < -M_PI / 2.0) && (phi1 > M_PI / 2.0)) {
            nphi2 = phi2 - 2.0 * M_PI;
        }

        Z = std::sqrt(std::abs(z1 * z2)) * std::exp(std::complex<double>(0.0, 1.0) * (nphi1 + nphi2) / 2.0);
    } else if (pasadas == 1) { // aritmetica
        Z = (z1 + z2) / 2.0;
    }

    return Z;
}

} // namespace farfield_m

