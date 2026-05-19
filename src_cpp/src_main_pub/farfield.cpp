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
