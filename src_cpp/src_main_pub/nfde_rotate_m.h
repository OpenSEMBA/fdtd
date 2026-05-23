#ifndef NFDE_ROTATE_M_H
#define NFDE_ROTATE_M_H

#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>

#include "nfde_types.h"

namespace nfde_rotate_m {

using NFDETypes_m::BloqueProbe_t;
using NFDETypes_m::coords_scaled_t;
using NFDETypes_m::coords_t;
using NFDETypes_m::Desplazamiento_t;
using NFDETypes_m::Dielectric_t;
using NFDETypes_m::FreqDepenMaterial_t;
using NFDETypes_m::MatrizMedios_t;
using NFDETypes_m::Parseador_t;
using NFDETypes_m::Sonda_t;

inline int sign_val(int val, int base) {
    return (base < 0) ? -val : val;
}

inline void rotateOrMpidir2(int32_t oor, int32_t& or_out) {
    if (oor == NFDETypes_m::iEx) or_out = NFDETypes_m::iEy;
    else if (oor == -NFDETypes_m::iEx) or_out = -NFDETypes_m::iEy;
    else if (oor == NFDETypes_m::iEy) or_out = NFDETypes_m::iEz;
    else if (oor == -NFDETypes_m::iEy) or_out = -NFDETypes_m::iEz;
    else if (oor == NFDETypes_m::iEz) or_out = NFDETypes_m::iEx;
    else if (oor == -NFDETypes_m::iEz) or_out = -NFDETypes_m::iEx;
}

inline void rotateOrMpidir1(int32_t oor, int32_t& or_out) {
    if (oor == NFDETypes_m::iEx) or_out = NFDETypes_m::iEz;
    else if (oor == -NFDETypes_m::iEx) or_out = -NFDETypes_m::iEz;
    else if (oor == NFDETypes_m::iEy) or_out = NFDETypes_m::iEx;
    else if (oor == -NFDETypes_m::iEy) or_out = -NFDETypes_m::iEx;
    else if (oor == NFDETypes_m::iEz) or_out = NFDETypes_m::iEy;
    else if (oor == -NFDETypes_m::iEz) or_out = -NFDETypes_m::iEy;
}

inline void ROTATEMPI(int mpidir, coords_t& coorden) {
    const int32_t oxi = coorden.Xi;
    const int32_t oxe = coorden.Xe;
    const int32_t oyi = coorden.Yi;
    const int32_t oye = coorden.Ye;
    const int32_t ozi = coorden.Zi;
    const int32_t oze = coorden.Ze;
    const int32_t txi = coorden.Xtrancos;
    const int32_t tyi = coorden.Ytrancos;
    const int32_t tzi = coorden.Ztrancos;
    const int32_t oor = coorden.Or;

    if (mpidir == 2) {
        coorden.Xi = ozi;
        coorden.Xe = oze;
        coorden.Xtrancos = tzi;
        coorden.Yi = oxi;
        coorden.Ye = oxe;
        coorden.Ytrancos = txi;
        coorden.Zi = oyi;
        coorden.Ze = oye;
        coorden.Ztrancos = tyi;
        rotateOrMpidir2(oor, coorden.Or);
    } else if (mpidir == 1) {
        coorden.Xi = oyi;
        coorden.Xe = oye;
        coorden.Xtrancos = tyi;
        coorden.Yi = ozi;
        coorden.Ye = oze;
        coorden.Ytrancos = tzi;
        coorden.Zi = oxi;
        coorden.Ze = oxe;
        coorden.Ztrancos = txi;
        rotateOrMpidir1(oor, coorden.Or);
    }
}

inline void ROTATEMPI_SCALED(int mpidir, coords_scaled_t& coorden) {
    const int32_t oxi = coorden.Xi;
    const int32_t oxe = coorden.Xe;
    const int32_t oyi = coorden.Yi;
    const int32_t oye = coorden.Ye;
    const int32_t ozi = coorden.Zi;
    const int32_t oze = coorden.Ze;
    const double oxc = coorden.xc;
    const double oyc = coorden.yc;
    const double ozc = coorden.zc;
    const int32_t oor = coorden.Or;

    if (mpidir == 2) {
        coorden.Xi = ozi;
        coorden.Xe = oze;
        coorden.Yi = oxi;
        coorden.Ye = oxe;
        coorden.Zi = oyi;
        coorden.Ze = oye;
        coorden.xc = ozc;
        coorden.yc = oxc;
        coorden.zc = oyc;
        rotateOrMpidir2(oor, coorden.Or);
    } else if (mpidir == 1) {
        coorden.Xi = oyi;
        coorden.Xe = oye;
        coorden.Yi = ozi;
        coorden.Ye = oze;
        coorden.Zi = oxi;
        coorden.Ze = oxe;
        coorden.xc = oyc;
        coorden.yc = ozc;
        coorden.zc = oxc;
        rotateOrMpidir1(oor, coorden.Or);
    }
}

inline void rotate_freq_depend_material_properties(int mpidir, FreqDepenMaterial_t& freqDepMat) {
    auto rotate_complex_block = [&](std::vector<std::complex<double>>& c11,
                                    std::vector<std::complex<double>>& c12,
                                    std::vector<std::complex<double>>& c13,
                                    std::vector<std::complex<double>>& c22,
                                    std::vector<std::complex<double>>& c23,
                                    std::vector<std::complex<double>>& c33) {
        if (c11.empty()) return;
        const std::complex<double> o11 = c11[0];
        const std::complex<double> o12 = c12[0];
        const std::complex<double> o13 = c13[0];
        const std::complex<double> o22 = c22[0];
        const std::complex<double> o23 = c23[0];
        const std::complex<double> o33 = c33[0];
        if (mpidir == 2) {
            c11[0] = o33;
            c12[0] = o23;
            c13[0] = o12;
            c22[0] = o11;
            c23[0] = o13;
            c33[0] = o22;
        } else if (mpidir == 1) {
            c11[0] = o22;
            c12[0] = o13;
            c13[0] = o23;
            c22[0] = o33;
            c23[0] = o12;
            c33[0] = o11;
        }
    };

    auto rotate_real_block = [&](double& r11, double& r12, double& r13,
                                 double& r22, double& r23, double& r33) {
        const double o11 = r11;
        const double o12 = r12;
        const double o13 = r13;
        const double o22 = r22;
        const double o23 = r23;
        const double o33 = r33;
        if (mpidir == 2) {
            r11 = o33;
            r12 = o23;
            r13 = o12;
            r22 = o11;
            r23 = o13;
            r33 = o22;
        } else if (mpidir == 1) {
            r11 = o22;
            r12 = o13;
            r13 = o23;
            r22 = o33;
            r23 = o12;
            r33 = o11;
        }
    };

    auto rotate_int_block = [&](int32_t& k11, int32_t& k12, int32_t& k13,
                                 int32_t& k22, int32_t& k23, int32_t& k33) {
        const int32_t o11 = k11;
        const int32_t o12 = k12;
        const int32_t o13 = k13;
        const int32_t o22 = k22;
        const int32_t o23 = k23;
        const int32_t o33 = k33;
        if (mpidir == 2) {
            k11 = o33;
            k12 = o23;
            k13 = o12;
            k22 = o11;
            k23 = o13;
            k33 = o22;
        } else if (mpidir == 1) {
            k11 = o22;
            k12 = o13;
            k13 = o23;
            k22 = o33;
            k23 = o12;
            k33 = o11;
        }
    };

    rotate_complex_block(freqDepMat.a11, freqDepMat.a12, freqDepMat.a13,
                         freqDepMat.a22, freqDepMat.a23, freqDepMat.a33);
    rotate_complex_block(freqDepMat.am11, freqDepMat.am12, freqDepMat.am13,
                         freqDepMat.am22, freqDepMat.am23, freqDepMat.am33);
    rotate_complex_block(freqDepMat.b11, freqDepMat.b12, freqDepMat.b13,
                         freqDepMat.b22, freqDepMat.b23, freqDepMat.b33);
    rotate_complex_block(freqDepMat.bm11, freqDepMat.bm12, freqDepMat.bm13,
                         freqDepMat.bm22, freqDepMat.bm23, freqDepMat.bm33);

    rotate_real_block(freqDepMat.eps11, freqDepMat.eps12, freqDepMat.eps13,
                      freqDepMat.eps22, freqDepMat.eps23, freqDepMat.eps33);
    rotate_real_block(freqDepMat.mu11, freqDepMat.mu12, freqDepMat.mu13,
                      freqDepMat.mu22, freqDepMat.mu23, freqDepMat.mu33);
    rotate_real_block(freqDepMat.sigma11, freqDepMat.sigma12, freqDepMat.sigma13,
                      freqDepMat.sigma22, freqDepMat.sigma23, freqDepMat.sigma33);
    rotate_real_block(freqDepMat.sigmam11, freqDepMat.sigmam12, freqDepMat.sigmam13,
                      freqDepMat.sigmam22, freqDepMat.sigmam23, freqDepMat.sigmam33);

    rotate_int_block(freqDepMat.K11, freqDepMat.K12, freqDepMat.K13,
                     freqDepMat.K22, freqDepMat.K23, freqDepMat.K33);
    rotate_int_block(freqDepMat.Km11, freqDepMat.Km12, freqDepMat.Km13,
                     freqDepMat.Km22, freqDepMat.Km23, freqDepMat.Km33);
}

inline void rotate_generateSpaceSteps(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.despl || !this_obj.matriz) return;

    const Desplazamiento_t old_despl = *this_obj.despl;
    const MatrizMedios_t old_matriz = *this_obj.matriz;

    const std::vector<double> poxi = old_despl.desX;
    const std::vector<double> poyi = old_despl.desY;
    const std::vector<double> pozi = old_despl.desZ;

    if (mpidir == 2) {
        this_obj.matriz->totalX = old_matriz.totalZ;
        this_obj.matriz->totalY = old_matriz.totalX;
        this_obj.matriz->totalZ = old_matriz.totalY;

        this_obj.despl->nX = old_despl.nZ;
        this_obj.despl->nY = old_despl.nX;
        this_obj.despl->nZ = old_despl.nY;

        this_obj.despl->mx1 = old_despl.mz1;
        this_obj.despl->my1 = old_despl.mx1;
        this_obj.despl->mz1 = old_despl.my1;

        this_obj.despl->mx2 = old_despl.mz2;
        this_obj.despl->my2 = old_despl.mx2;
        this_obj.despl->mz2 = old_despl.my2;

        this_obj.despl->originx = old_despl.originz;
        this_obj.despl->originy = old_despl.originx;
        this_obj.despl->originz = old_despl.originy;

        this_obj.despl->desX = pozi;
        this_obj.despl->desY = poxi;
        this_obj.despl->desZ = poyi;
    } else if (mpidir == 1) {
        this_obj.matriz->totalX = old_matriz.totalY;
        this_obj.matriz->totalY = old_matriz.totalZ;
        this_obj.matriz->totalZ = old_matriz.totalX;

        this_obj.despl->nX = old_despl.nY;
        this_obj.despl->nY = old_despl.nZ;
        this_obj.despl->nZ = old_despl.nX;

        this_obj.despl->mx1 = old_despl.my1;
        this_obj.despl->my1 = old_despl.mz1;
        this_obj.despl->mz1 = old_despl.mx1;

        this_obj.despl->mx2 = old_despl.my2;
        this_obj.despl->my2 = old_despl.mz2;
        this_obj.despl->mz2 = old_despl.mx2;

        this_obj.despl->originx = old_despl.originy;
        this_obj.despl->originy = old_despl.originz;
        this_obj.despl->originz = old_despl.originx;

        this_obj.despl->desX = poyi;
        this_obj.despl->desY = pozi;
        this_obj.despl->desZ = poxi;
    }
}

inline void rotate_generateCurrent_Field_Sources(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.nodSrc) return;
    const int tama = this_obj.nodSrc->n_nodSrc;
    for (int i = 0; i < tama; ++i) {
        const int tama2 = this_obj.nodSrc->NodalSource[i].n_C1P;
        for (int ii = 0; ii < tama2; ++ii) {
            ROTATEMPI_SCALED(mpidir, this_obj.nodSrc->NodalSource[i].c1P[ii]);
        }
        const int tama3 = this_obj.nodSrc->NodalSource[i].n_C2P;
        for (int ii = 0; ii < tama3; ++ii) {
            ROTATEMPI_SCALED(mpidir, this_obj.nodSrc->NodalSource[i].c2P[ii]);
        }
    }
}

inline void rotate_generatePlaneWaves(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.plnSrc) return;
    const int tama = this_obj.plnSrc->nc;
    for (int i = 0; i < tama; ++i) {
        const double theta = this_obj.plnSrc->collection[i].theta;
        const double phi = this_obj.plnSrc->collection[i].phi;
        const double alpha = this_obj.plnSrc->collection[i].alpha;
        const double beta = this_obj.plnSrc->collection[i].beta;

        if (mpidir == 2) {
            const int oxi = this_obj.plnSrc->collection[i].coor1[0];
            const int oxe = this_obj.plnSrc->collection[i].coor2[0];
            const int oyi = this_obj.plnSrc->collection[i].coor1[1];
            const int oye = this_obj.plnSrc->collection[i].coor2[1];
            const int ozi = this_obj.plnSrc->collection[i].coor1[2];
            const int oze = this_obj.plnSrc->collection[i].coor2[2];

            this_obj.plnSrc->collection[i].coor1[0] = ozi;
            this_obj.plnSrc->collection[i].coor2[0] = oze;
            this_obj.plnSrc->collection[i].coor1[1] = oxi;
            this_obj.plnSrc->collection[i].coor2[1] = oxe;
            this_obj.plnSrc->collection[i].coor1[2] = oyi;
            this_obj.plnSrc->collection[i].coor2[2] = oye;

            this_obj.plnSrc->collection[i].theta = std::atan2(
                std::sqrt(std::cos(theta) * std::cos(theta) +
                          std::cos(phi) * std::cos(phi) * std::sin(theta) * std::sin(theta)),
                std::sin(phi) * std::sin(theta));
            this_obj.plnSrc->collection[i].phi =
                std::atan2(std::cos(phi) * std::sin(theta), std::cos(theta));
            this_obj.plnSrc->collection[i].alpha = std::atan2(
                std::sqrt(std::cos(alpha) * std::cos(alpha) +
                          std::cos(beta) * std::cos(beta) * std::sin(alpha) * std::sin(alpha)),
                std::sin(beta) * std::sin(alpha));
            this_obj.plnSrc->collection[i].beta =
                std::atan2(std::cos(beta) * std::sin(alpha), std::cos(alpha));
        } else if (mpidir == 1) {
            const int oxi = this_obj.plnSrc->collection[i].coor1[0];
            const int oxe = this_obj.plnSrc->collection[i].coor2[0];
            const int oyi = this_obj.plnSrc->collection[i].coor1[1];
            const int oye = this_obj.plnSrc->collection[i].coor2[1];
            const int ozi = this_obj.plnSrc->collection[i].coor1[2];
            const int oze = this_obj.plnSrc->collection[i].coor2[2];

            this_obj.plnSrc->collection[i].coor1[0] = oyi;
            this_obj.plnSrc->collection[i].coor2[0] = oye;
            this_obj.plnSrc->collection[i].coor1[1] = ozi;
            this_obj.plnSrc->collection[i].coor2[1] = oze;
            this_obj.plnSrc->collection[i].coor1[2] = oxi;
            this_obj.plnSrc->collection[i].coor2[2] = oxe;

            this_obj.plnSrc->collection[i].theta = std::atan2(
                std::sqrt(std::cos(theta) * std::cos(theta) +
                          std::sin(phi) * std::sin(phi) * std::sin(theta) * std::sin(theta)),
                std::cos(phi) * std::sin(theta));
            this_obj.plnSrc->collection[i].phi =
                std::atan2(std::cos(theta), std::sin(phi) * std::sin(theta));
            this_obj.plnSrc->collection[i].alpha = std::atan2(
                std::sqrt(std::cos(alpha) * std::cos(alpha) +
                          std::sin(beta) * std::sin(beta) * std::sin(alpha) * std::sin(alpha)),
                std::cos(beta) * std::sin(alpha));
            this_obj.plnSrc->collection[i].beta =
                std::atan2(std::cos(alpha), std::sin(beta) * std::sin(alpha));
        }
    }
}

inline void rotate_generateBoxSources(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.boxSrc) return;
    const int tama = this_obj.boxSrc->nVols;
    for (int i = 0; i < tama; ++i) {
        const int oxi = this_obj.boxSrc->Vols[i].coor1[0];
        const int oxe = this_obj.boxSrc->Vols[i].coor2[0];
        const int oyi = this_obj.boxSrc->Vols[i].coor1[1];
        const int oye = this_obj.boxSrc->Vols[i].coor2[1];
        const int ozi = this_obj.boxSrc->Vols[i].coor1[2];
        const int oze = this_obj.boxSrc->Vols[i].coor2[2];

        if (mpidir == 2) {
            this_obj.boxSrc->Vols[i].coor1[0] = ozi;
            this_obj.boxSrc->Vols[i].coor2[0] = oze;
            this_obj.boxSrc->Vols[i].coor1[1] = oxi;
            this_obj.boxSrc->Vols[i].coor2[1] = oxe;
            this_obj.boxSrc->Vols[i].coor1[2] = oyi;
            this_obj.boxSrc->Vols[i].coor2[2] = oye;
        } else if (mpidir == 1) {
            this_obj.boxSrc->Vols[i].coor1[0] = oyi;
            this_obj.boxSrc->Vols[i].coor2[0] = oye;
            this_obj.boxSrc->Vols[i].coor1[1] = ozi;
            this_obj.boxSrc->Vols[i].coor2[1] = oze;
            this_obj.boxSrc->Vols[i].coor1[2] = oxi;
            this_obj.boxSrc->Vols[i].coor2[2] = oxe;
        }
    }
}

inline void rotate_generateFronteras(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.front) return;

    const int oxl = this_obj.front->tipoFrontera[0];
    const int oxu = this_obj.front->tipoFrontera[1];
    const int oyl = this_obj.front->tipoFrontera[2];
    const int oyu = this_obj.front->tipoFrontera[3];
    const int ozl = this_obj.front->tipoFrontera[4];
    const int ozu = this_obj.front->tipoFrontera[5];

    if (mpidir == 2) {
        this_obj.front->tipoFrontera[0] = ozl;
        this_obj.front->tipoFrontera[1] = ozu;
        this_obj.front->tipoFrontera[2] = oxl;
        this_obj.front->tipoFrontera[3] = oxu;
        this_obj.front->tipoFrontera[4] = oyl;
        this_obj.front->tipoFrontera[5] = oyu;
    } else if (mpidir == 1) {
        this_obj.front->tipoFrontera[0] = oyl;
        this_obj.front->tipoFrontera[1] = oyu;
        this_obj.front->tipoFrontera[2] = ozl;
        this_obj.front->tipoFrontera[3] = ozu;
        this_obj.front->tipoFrontera[4] = oxl;
        this_obj.front->tipoFrontera[5] = oxu;
    }

    if (mpidir == 1 || mpidir == 2) {
        const double orden_xl = this_obj.front->propiedadesPML[0].orden;
        const double orden_xu = this_obj.front->propiedadesPML[1].orden;
        const double orden_yl = this_obj.front->propiedadesPML[2].orden;
        const double orden_yu = this_obj.front->propiedadesPML[3].orden;
        const double orden_zl = this_obj.front->propiedadesPML[4].orden;
        const double orden_zu = this_obj.front->propiedadesPML[5].orden;

        const double refl_xl = this_obj.front->propiedadesPML[0].refl;
        const double refl_xu = this_obj.front->propiedadesPML[1].refl;
        const double refl_yl = this_obj.front->propiedadesPML[2].refl;
        const double refl_yu = this_obj.front->propiedadesPML[3].refl;
        const double refl_zl = this_obj.front->propiedadesPML[4].refl;
        const double refl_zu = this_obj.front->propiedadesPML[5].refl;

        const int32_t numcapas_xl = this_obj.front->propiedadesPML[0].numCapas;
        const int32_t numcapas_xu = this_obj.front->propiedadesPML[1].numCapas;
        const int32_t numcapas_yl = this_obj.front->propiedadesPML[2].numCapas;
        const int32_t numcapas_yu = this_obj.front->propiedadesPML[3].numCapas;
        const int32_t numcapas_zl = this_obj.front->propiedadesPML[4].numCapas;
        const int32_t numcapas_zu = this_obj.front->propiedadesPML[5].numCapas;

        if (mpidir == 2) {
            this_obj.front->propiedadesPML[0].orden = orden_zl;
            this_obj.front->propiedadesPML[1].orden = orden_zu;
            this_obj.front->propiedadesPML[2].orden = orden_xl;
            this_obj.front->propiedadesPML[3].orden = orden_xu;
            this_obj.front->propiedadesPML[4].orden = orden_yl;
            this_obj.front->propiedadesPML[5].orden = orden_yu;

            this_obj.front->propiedadesPML[0].refl = refl_zl;
            this_obj.front->propiedadesPML[1].refl = refl_zu;
            this_obj.front->propiedadesPML[2].refl = refl_xl;
            this_obj.front->propiedadesPML[3].refl = refl_xu;
            this_obj.front->propiedadesPML[4].refl = refl_yl;
            this_obj.front->propiedadesPML[5].refl = refl_yu;

            this_obj.front->propiedadesPML[0].numCapas = numcapas_zl;
            this_obj.front->propiedadesPML[1].numCapas = numcapas_zu;
            this_obj.front->propiedadesPML[2].numCapas = numcapas_xl;
            this_obj.front->propiedadesPML[3].numCapas = numcapas_xu;
            this_obj.front->propiedadesPML[4].numCapas = numcapas_yl;
            this_obj.front->propiedadesPML[5].numCapas = numcapas_yu;
        } else {
            this_obj.front->propiedadesPML[0].orden = orden_yl;
            this_obj.front->propiedadesPML[1].orden = orden_yu;
            this_obj.front->propiedadesPML[2].orden = orden_zl;
            this_obj.front->propiedadesPML[3].orden = orden_zu;
            this_obj.front->propiedadesPML[4].orden = orden_xl;
            this_obj.front->propiedadesPML[5].orden = orden_xu;

            this_obj.front->propiedadesPML[0].refl = refl_yl;
            this_obj.front->propiedadesPML[1].refl = refl_yu;
            this_obj.front->propiedadesPML[2].refl = refl_zl;
            this_obj.front->propiedadesPML[3].refl = refl_zu;
            this_obj.front->propiedadesPML[4].refl = refl_xl;
            this_obj.front->propiedadesPML[5].refl = refl_xu;

            this_obj.front->propiedadesPML[0].numCapas = numcapas_yl;
            this_obj.front->propiedadesPML[1].numCapas = numcapas_yu;
            this_obj.front->propiedadesPML[2].numCapas = numcapas_zl;
            this_obj.front->propiedadesPML[3].numCapas = numcapas_zu;
            this_obj.front->propiedadesPML[4].numCapas = numcapas_xl;
            this_obj.front->propiedadesPML[5].numCapas = numcapas_xu;
        }
    }
}

inline void rotate_generatePECs(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.pecRegs) return;
    for (int i = 0; i < this_obj.pecRegs->nVols; ++i) {
        ROTATEMPI(mpidir, this_obj.pecRegs->Vols[i]);
    }
    for (int i = 0; i < this_obj.pecRegs->nSurfs; ++i) {
        ROTATEMPI(mpidir, this_obj.pecRegs->Surfs[i]);
    }
    for (int i = 0; i < this_obj.pecRegs->nLins; ++i) {
        ROTATEMPI(mpidir, this_obj.pecRegs->Lins[i]);
    }
}

inline void rotate_dielectric(Dielectric_t& body, int mpidir) {
    for (int ii = 0; ii < body.n_C1P; ++ii) {
        ROTATEMPI(mpidir, body.c1P[ii]);
    }
    if (body.n_C1P > 0) {
        body.DiodOri = body.c1P[body.n_C1P - 1].Or;
    }
    for (int ii = 0; ii < body.n_C2P; ++ii) {
        ROTATEMPI(mpidir, body.c2P[ii]);
    }
    if (body.n_C2P > 0) {
        body.DiodOri = body.c2P[body.n_C2P - 1].Or;
    }
}

inline void rotate_generateNONMetals(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.DielRegs) return;
    for (int i = 0; i < this_obj.DielRegs->nVols; ++i) {
        rotate_dielectric(this_obj.DielRegs->Vols[i], mpidir);
    }
    for (int i = 0; i < this_obj.DielRegs->nSurfs; ++i) {
        rotate_dielectric(this_obj.DielRegs->Surfs[i], mpidir);
    }
    for (int i = 0; i < this_obj.DielRegs->nLins; ++i) {
        rotate_dielectric(this_obj.DielRegs->Lins[i], mpidir);
    }
}

inline void rotate_thin_wire_dir(int mpidir, int32_t d, int32_t& out_d) {
    if (mpidir == 2) {
        if (d == NFDETypes_m::iEx) out_d = NFDETypes_m::iEy;
        else if (d == NFDETypes_m::iEy) out_d = NFDETypes_m::iEz;
        else if (d == NFDETypes_m::iEz) out_d = NFDETypes_m::iEx;
        else out_d = d;
    } else if (mpidir == 1) {
        if (d == NFDETypes_m::iEx) out_d = NFDETypes_m::iEz;
        else if (d == NFDETypes_m::iEy) out_d = NFDETypes_m::iEx;
        else if (d == NFDETypes_m::iEz) out_d = NFDETypes_m::iEy;
        else out_d = d;
    } else {
        out_d = d;
    }
}

inline void rotate_generateThinWires(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.tWires) return;
    const int tama = this_obj.tWires->n_tw;
    for (int i = 0; i < tama; ++i) {
        const int tama2 = this_obj.tWires->tw[i].n_twc;
        for (int ii = 0; ii < tama2; ++ii) {
            const int oldx = this_obj.tWires->tw[i].twc[ii].i;
            const int oldy = this_obj.tWires->tw[i].twc[ii].j;
            const int oldz = this_obj.tWires->tw[i].twc[ii].K;
            const int32_t old_d = this_obj.tWires->tw[i].twc[ii].d;

            if (mpidir == 2) {
                this_obj.tWires->tw[i].twc[ii].i = oldz;
                this_obj.tWires->tw[i].twc[ii].j = oldx;
                this_obj.tWires->tw[i].twc[ii].K = oldy;
            } else if (mpidir == 1) {
                this_obj.tWires->tw[i].twc[ii].i = oldy;
                this_obj.tWires->tw[i].twc[ii].j = oldz;
                this_obj.tWires->tw[i].twc[ii].K = oldx;
            }
            rotate_thin_wire_dir(mpidir, old_d, this_obj.tWires->tw[i].twc[ii].d);
        }
    }
}

inline void rotate_generateSlantedWires(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.sWires) return;
    const int tama = this_obj.sWires->n_sw;
    for (int i = 0; i < tama; ++i) {
        const int tama2 = this_obj.sWires->sw[i].n_swc;
        for (int ii = 0; ii < tama2; ++ii) {
            const double oldx = this_obj.sWires->sw[i].swc[ii].x;
            const double oldy = this_obj.sWires->sw[i].swc[ii].y;
            const double oldz = this_obj.sWires->sw[i].swc[ii].z;

            if (mpidir == 2) {
                this_obj.sWires->sw[i].swc[ii].x = oldz;
                this_obj.sWires->sw[i].swc[ii].y = oldx;
                this_obj.sWires->sw[i].swc[ii].z = oldy;
            } else if (mpidir == 1) {
                this_obj.sWires->sw[i].swc[ii].x = oldy;
                this_obj.sWires->sw[i].swc[ii].y = oldz;
                this_obj.sWires->sw[i].swc[ii].z = oldx;
            }
        }
    }
}

inline void rotate_generateThinSlots(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.tSlots) return;
    const int tama = this_obj.tSlots->n_tg;
    for (int i = 0; i < tama; ++i) {
        const int tama2 = this_obj.tSlots->tg[i].n_tgc;
        for (int ii = 0; ii < tama2; ++ii) {
            const int oldx = this_obj.tSlots->tg[i].tgc[ii].i;
            const int oldy = this_obj.tSlots->tg[i].tgc[ii].j;
            const int oldz = this_obj.tSlots->tg[i].tgc[ii].K;
            const int32_t old_dir = this_obj.tSlots->tg[i].tgc[ii].dir;

            if (mpidir == 2) {
                this_obj.tSlots->tg[i].tgc[ii].i = oldz;
                this_obj.tSlots->tg[i].tgc[ii].j = oldx;
                this_obj.tSlots->tg[i].tgc[ii].K = oldy;
            } else if (mpidir == 1) {
                this_obj.tSlots->tg[i].tgc[ii].i = oldy;
                this_obj.tSlots->tg[i].tgc[ii].j = oldz;
                this_obj.tSlots->tg[i].tgc[ii].K = oldx;
            }
            rotate_thin_wire_dir(mpidir, old_dir, this_obj.tSlots->tg[i].tgc[ii].dir);
        }
    }
}

inline void rotate_generateLossyThinSurface(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.LossyThinSurfs) return;
    const int tama = this_obj.LossyThinSurfs->length;
    for (int i = 0; i < tama; ++i) {
        const int tama2 = this_obj.LossyThinSurfs->cs[i].nc;
        for (int ii = 0; ii < tama2; ++ii) {
            ROTATEMPI(mpidir, this_obj.LossyThinSurfs->cs[i].c[ii]);
        }
    }
}

inline void rotate_generateFDMs(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.frqDepMats) return;

    auto rotate_fdm_body = [&](FreqDepenMaterial_t& body) {
        for (int ii = 0; ii < body.n_c; ++ii) {
            ROTATEMPI(mpidir, body.c[ii]);
        }
        rotate_freq_depend_material_properties(mpidir, body);
    };

    for (int i = 0; i < this_obj.frqDepMats->nVols; ++i) {
        rotate_fdm_body(this_obj.frqDepMats->Vols[i]);
    }
    for (int i = 0; i < this_obj.frqDepMats->nSurfs; ++i) {
        rotate_fdm_body(this_obj.frqDepMats->Surfs[i]);
    }
    for (int i = 0; i < this_obj.frqDepMats->nLins; ++i) {
        rotate_fdm_body(this_obj.frqDepMats->Lins[i]);
    }
}

inline void rotate_sonda_coords(int mpidir, Sonda_t& probe) {
    for (int iii = 0; iii < probe.n_cord; ++iii) {
        if (mpidir == 2) {
            const int iox = probe.i[iii];
            const int ioy = probe.j[iii];
            const int ioz = probe.K[iii];
            probe.i[iii] = ioz;
            probe.j[iii] = iox;
            probe.K[iii] = ioy;
        } else if (mpidir == 1) {
            const int iox = probe.i[iii];
            const int ioy = probe.j[iii];
            const int ioz = probe.K[iii];
            probe.i[iii] = ioy;
            probe.j[iii] = ioz;
            probe.K[iii] = iox;
        }
    }
}

inline void rotate_generateSONDAs(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.oldSONDA) return;
    const int tama = this_obj.oldSONDA->n_probes;
    for (int i = 0; i < tama; ++i) {
        for (int ii = 0; ii < this_obj.oldSONDA->probes[i].n_FarField; ++ii) {
            const double thetastart = this_obj.oldSONDA->probes[i].FarField[ii].probe.thetastart;
            const double thetastop = this_obj.oldSONDA->probes[i].FarField[ii].probe.thetastop;
            const double phistart = this_obj.oldSONDA->probes[i].FarField[ii].probe.phistart;
            const double phistop = this_obj.oldSONDA->probes[i].FarField[ii].probe.phistop;

            if (mpidir == 2) {
                this_obj.oldSONDA->probes[i].FarField[ii].probe.thetastart = std::atan2(
                    std::sqrt(std::cos(thetastart) * std::cos(thetastart) +
                              std::cos(phistart) * std::cos(phistart) * std::sin(thetastart) *
                                  std::sin(thetastart)),
                    std::sin(phistart) * std::sin(thetastart));
                this_obj.oldSONDA->probes[i].FarField[ii].probe.phistart =
                    std::atan2(std::cos(phistart) * std::sin(thetastart), std::cos(thetastart));
                this_obj.oldSONDA->probes[i].FarField[ii].probe.thetastop = std::atan2(
                    std::sqrt(std::cos(thetastop) * std::cos(thetastop) +
                              std::cos(phistop) * std::cos(phistop) * std::sin(thetastop) *
                                  std::sin(thetastop)),
                    std::sin(phistop) * std::sin(thetastop));
                this_obj.oldSONDA->probes[i].FarField[ii].probe.phistop =
                    std::atan2(std::cos(phistop) * std::sin(thetastop), std::cos(thetastop));
            } else if (mpidir == 1) {
                this_obj.oldSONDA->probes[i].FarField[ii].probe.thetastart = std::atan2(
                    std::sqrt(std::cos(thetastart) * std::cos(thetastart) +
                              std::sin(phistart) * std::sin(phistart) * std::sin(thetastart) *
                                  std::sin(thetastart)),
                    std::cos(phistart) * std::sin(thetastart));
                this_obj.oldSONDA->probes[i].FarField[ii].probe.phistart =
                    std::atan2(std::cos(thetastart), std::sin(phistart) * std::sin(thetastart));
                this_obj.oldSONDA->probes[i].FarField[ii].probe.thetastop = std::atan2(
                    std::sqrt(std::cos(thetastop) * std::cos(thetastop) +
                              std::sin(phistop) * std::sin(phistop) * std::sin(thetastop) *
                                  std::sin(thetastop)),
                    std::cos(phistop) * std::sin(thetastop));
                this_obj.oldSONDA->probes[i].FarField[ii].probe.phistop =
                    std::atan2(std::cos(thetastop), std::sin(phistop) * std::sin(thetastop));
            }
            rotate_sonda_coords(mpidir, this_obj.oldSONDA->probes[i].FarField[ii].probe);
        }
    }

    for (int i = 0; i < tama; ++i) {
        for (int ii = 0; ii < this_obj.oldSONDA->probes[i].n_Electric; ++ii) {
            rotate_sonda_coords(mpidir, this_obj.oldSONDA->probes[i].Electric[ii].probe);
        }
    }
    for (int i = 0; i < tama; ++i) {
        for (int ii = 0; ii < this_obj.oldSONDA->probes[i].n_Magnetic; ++ii) {
            rotate_sonda_coords(mpidir, this_obj.oldSONDA->probes[i].Magnetic[ii].probe);
        }
    }
}

inline bool is_massonda_rotatable(int32_t oor) {
    return oor == NFDETypes_m::NP_COR_EX || oor == NFDETypes_m::NP_COR_EY ||
           oor == NFDETypes_m::NP_COR_EZ || oor == NFDETypes_m::NP_COR_HX ||
           oor == NFDETypes_m::NP_COR_HY || oor == NFDETypes_m::NP_COR_HZ;
}

inline void rotate_massonda_or(int mpidir, int32_t oor, int32_t& or_out) {
    if (mpidir == 2) {
        if (oor == NFDETypes_m::NP_COR_EX) or_out = NFDETypes_m::NP_COR_EY;
        else if (oor == NFDETypes_m::NP_COR_EY) or_out = NFDETypes_m::NP_COR_EZ;
        else if (oor == NFDETypes_m::NP_COR_EZ) or_out = NFDETypes_m::NP_COR_EX;
        else if (oor == NFDETypes_m::NP_COR_HX) or_out = NFDETypes_m::NP_COR_HY;
        else if (oor == NFDETypes_m::NP_COR_HY) or_out = NFDETypes_m::NP_COR_HZ;
        else if (oor == NFDETypes_m::NP_COR_HZ) or_out = NFDETypes_m::NP_COR_HX;
    } else if (mpidir == 1) {
        if (oor == NFDETypes_m::NP_COR_EX) or_out = NFDETypes_m::NP_COR_EZ;
        else if (oor == NFDETypes_m::NP_COR_EY) or_out = NFDETypes_m::NP_COR_EX;
        else if (oor == NFDETypes_m::NP_COR_EZ) or_out = NFDETypes_m::NP_COR_EY;
        else if (oor == NFDETypes_m::NP_COR_HX) or_out = NFDETypes_m::NP_COR_HZ;
        else if (oor == NFDETypes_m::NP_COR_HY) or_out = NFDETypes_m::NP_COR_HX;
        else if (oor == NFDETypes_m::NP_COR_HZ) or_out = NFDETypes_m::NP_COR_HY;
    }
}

inline void rotate_generateMasSondas(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.Sonda) return;
    const int tama = this_obj.Sonda->length;
    for (int i = 0; i < tama; ++i) {
        const int tama2 = this_obj.Sonda->collection[i].len_cor;
        for (int ii = 0; ii < tama2; ++ii) {
            coords_t& coord = this_obj.Sonda->collection[i].cordinates[ii];
            const int32_t oxi = coord.Xi;
            const int32_t oxe = coord.Xe;
            const int32_t oyi = coord.Yi;
            const int32_t oye = coord.Ye;
            const int32_t ozi = coord.Zi;
            const int32_t oze = coord.Ze;
            const int32_t oor = coord.Or;
            const int32_t txi = coord.Xtrancos;
            const int32_t tyi = coord.Ytrancos;
            const int32_t tzi = coord.Ztrancos;

            if (!is_massonda_rotatable(oor)) return;

            if (mpidir == 2) {
                coord.Xi = ozi;
                coord.Xe = oze;
                coord.Xtrancos = tzi;
                coord.Yi = oxi;
                coord.Ye = oxe;
                coord.Ytrancos = txi;
                coord.Zi = oyi;
                coord.Ze = oye;
                coord.Ztrancos = tyi;
            } else if (mpidir == 1) {
                coord.Xi = oyi;
                coord.Xe = oye;
                coord.Xtrancos = tyi;
                coord.Yi = ozi;
                coord.Ye = oze;
                coord.Ytrancos = tzi;
                coord.Zi = oxi;
                coord.Ze = oxe;
                coord.Ztrancos = txi;
            }
            rotate_massonda_or(mpidir, oor, coord.Or);
        }
    }
}

inline void rotate_generateBloqueProbes(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.BloquePrb) return;
    const int tama = this_obj.BloquePrb->n_bp;
    for (int i = 0; i < tama; ++i) {
        BloqueProbe_t& bp = this_obj.BloquePrb->bp[i];
        const int oxi = bp.i1;
        const int oxe = bp.i2;
        const int oyi = bp.j1;
        const int oye = bp.j2;
        const int ozi = bp.k1;
        const int oze = bp.k2;
        const int32_t nml = bp.nml;

        if (mpidir == 2) {
            bp.i1 = ozi;
            bp.i2 = oze;
            bp.j1 = oxi;
            bp.j2 = oxe;
            bp.k1 = oyi;
            bp.k2 = oye;
            rotate_thin_wire_dir(mpidir, nml, bp.nml);
        } else if (mpidir == 1) {
            bp.i1 = oyi;
            bp.i2 = oye;
            bp.j1 = ozi;
            bp.j2 = oze;
            bp.k1 = oxi;
            bp.k2 = oxe;
            rotate_thin_wire_dir(mpidir, nml, bp.nml);
        }
    }
}

inline bool is_volumic_rotatable(int32_t oor) {
    return oor == NFDETypes_m::iExC || oor == NFDETypes_m::iEyC || oor == NFDETypes_m::iEzC ||
           oor == NFDETypes_m::iHxC || oor == NFDETypes_m::iHyC || oor == NFDETypes_m::iHzC ||
           oor == NFDETypes_m::iCurX || oor == NFDETypes_m::iCurY || oor == NFDETypes_m::iCurZ ||
           oor == NFDETypes_m::iMEC || oor == NFDETypes_m::iMHC || oor == NFDETypes_m::iCur;
}

inline void rotate_volumic_or(int mpidir, int32_t oor, int32_t& or_out) {
    if (mpidir == 2) {
        if (oor == NFDETypes_m::iExC) or_out = NFDETypes_m::iEyC;
        else if (oor == NFDETypes_m::iEyC) or_out = NFDETypes_m::iEzC;
        else if (oor == NFDETypes_m::iEzC) or_out = NFDETypes_m::iExC;
        else if (oor == NFDETypes_m::iHxC) or_out = NFDETypes_m::iHyC;
        else if (oor == NFDETypes_m::iHyC) or_out = NFDETypes_m::iHzC;
        else if (oor == NFDETypes_m::iHzC) or_out = NFDETypes_m::iHxC;
        else if (oor == NFDETypes_m::iCurX) or_out = NFDETypes_m::iCurY;
        else if (oor == NFDETypes_m::iCurY) or_out = NFDETypes_m::iCurZ;
        else if (oor == NFDETypes_m::iCurZ) or_out = NFDETypes_m::iCurX;
    } else if (mpidir == 1) {
        if (oor == NFDETypes_m::iExC) or_out = NFDETypes_m::iEzC;
        else if (oor == NFDETypes_m::iEyC) or_out = NFDETypes_m::iExC;
        else if (oor == NFDETypes_m::iEzC) or_out = NFDETypes_m::iEyC;
        else if (oor == NFDETypes_m::iHxC) or_out = NFDETypes_m::iHzC;
        else if (oor == NFDETypes_m::iHyC) or_out = NFDETypes_m::iHxC;
        else if (oor == NFDETypes_m::iHzC) or_out = NFDETypes_m::iHyC;
        else if (oor == NFDETypes_m::iCurX) or_out = NFDETypes_m::iCurZ;
        else if (oor == NFDETypes_m::iCurY) or_out = NFDETypes_m::iCurX;
        else if (oor == NFDETypes_m::iCurZ) or_out = NFDETypes_m::iCurY;
    }
}

inline void rotate_generateVolumicProbes(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.VolPrb) return;
    const int tama = this_obj.VolPrb->length;
    for (int i = 0; i < tama; ++i) {
        const int tama2 = this_obj.VolPrb->collection[i].len_cor;
        for (int ii = 0; ii < tama2; ++ii) {
            coords_t& coord = this_obj.VolPrb->collection[i].cordinates[ii];
            const int32_t oxi = coord.Xi;
            const int32_t oxe = coord.Xe;
            const int32_t oyi = coord.Yi;
            const int32_t oye = coord.Ye;
            const int32_t ozi = coord.Zi;
            const int32_t oze = coord.Ze;
            const int32_t oor = coord.Or;
            const int32_t txi = coord.Xtrancos;
            const int32_t tyi = coord.Ytrancos;
            const int32_t tzi = coord.Ztrancos;

            if (!is_volumic_rotatable(oor)) return;

            if (mpidir == 2) {
                coord.Xi = ozi;
                coord.Xe = oze;
                coord.Yi = oxi;
                coord.Ye = oxe;
                coord.Zi = oyi;
                coord.Ze = oye;
                coord.Xtrancos = tzi;
                coord.Ytrancos = txi;
                coord.Ztrancos = tyi;
            } else if (mpidir == 1) {
                coord.Xi = oyi;
                coord.Xe = oye;
                coord.Yi = ozi;
                coord.Ye = oze;
                coord.Zi = oxi;
                coord.Ze = oxe;
                coord.Xtrancos = tyi;
                coord.Ytrancos = tzi;
                coord.Ztrancos = txi;
            }
            rotate_volumic_or(mpidir, oor, coord.Or);
        }
    }
}

#ifdef CompileWithMTLN
inline void rotate_mtln(Parseador_t& this_obj, int mpidir) {
    if (!this_obj.mtln) return;
    for (size_t i = 0; i < this_obj.mtln->cables.size(); ++i) {
        auto* cable = this_obj.mtln->cables[i].ptr.get();
        if (!cable) continue;
        for (size_t j = 0; j < cable->segments.size(); ++j) {
            const int x = cable->segments[j].x;
            const int y = cable->segments[j].y;
            const int z = cable->segments[j].z;
            const int or_val = cable->segments[j].orientation;

            if (mpidir == 2) {
                cable->segments[j].x = z;
                cable->segments[j].y = x;
                cable->segments[j].z = y;
                switch (std::abs(or_val)) {
                    case 1: cable->segments[j].orientation = sign_val(2, or_val); break;
                    case 2: cable->segments[j].orientation = sign_val(3, or_val); break;
                    case 3: cable->segments[j].orientation = sign_val(1, or_val); break;
                    default: break;
                }
            } else if (mpidir == 1) {
                cable->segments[j].x = y;
                cable->segments[j].y = z;
                cable->segments[j].z = x;
                switch (std::abs(or_val)) {
                    case 1: cable->segments[j].orientation = sign_val(3, or_val); break;
                    case 2: cable->segments[j].orientation = sign_val(1, or_val); break;
                    case 3: cable->segments[j].orientation = sign_val(2, or_val); break;
                    default: break;
                }
            }
        }
    }
}
#endif

} // namespace nfde_rotate_m

#endif
