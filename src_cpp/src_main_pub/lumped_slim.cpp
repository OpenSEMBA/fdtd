#include "lumped_slim.h"

#include <cmath>
#include <iostream>

namespace lumped_slim {

void LumpedSolver_t::calcConstants(double dt, double eps0, double mu0) {
    (void)mu0;
    for (auto& lumped : nodes) {
        const auto& m = lumped.mat;
        const double epsilon = m.epr;
        const double sigma = m.sigma;
        const double Resist = m.R;
        const double Induct = m.L;
        const double Capaci = m.C;

        const double alignedDeltaE = lumped.alignedDeltaE;
        const double transversalDeltaHa = lumped.transversalDeltaHa;
        const double transversalDeltaHb = lumped.transversalDeltaHb;

        const double epsilonEffCapac =
            alignedDeltaE * Capaci / (transversalDeltaHa * transversalDeltaHb);
        const double sigmaEffResistInduct =
            alignedDeltaE * dt /
            (2.0 * transversalDeltaHa * transversalDeltaHb * (Induct + Resist * dt / 2.0));
        const double sigmaEffResist =
            alignedDeltaE / (Resist * transversalDeltaHa * transversalDeltaHb);
        const double sigmaEffResistCapac = sigmaEffResist;
        const double sigmaEffResistDiode = sigmaEffResist;
        const double currentCoeff =
            (Induct - Resist * dt / 2.0) / (Induct + Resist * dt / 2.0);

        double sigmaeff = 0.0;
        double epsiloneff = 0.0;
        if (m.resistor) {
            sigmaeff = sigma + sigmaEffResist;
            epsiloneff = epsilon;
        } else if (m.inductor) {
            sigmaeff = sigma + sigmaEffResistInduct;
            epsiloneff = epsilon;
        } else if (m.capacitor) {
            sigmaeff = sigma + sigmaEffResistCapac;
            epsiloneff = epsilon + epsilonEffCapac;
        } else if (m.diodo) {
            sigmaeff = sigma + sigmaEffResistDiode;
            epsiloneff = epsilon;
        }

        const double G1 =
            (epsiloneff / dt - sigmaeff / 2.0) / (epsiloneff / dt + sigmaeff / 2.0);
        const double G2 = 1.0 / (epsiloneff / dt + sigmaeff / 2.0);

        lumped.G1 = G1;
        lumped.G2a = G2 / transversalDeltaHa;
        lumped.G2b = G2 / transversalDeltaHb;
        lumped.GJ = G2 * (1.0 + currentCoeff) / 2.0;
        lumped.sigmaEffResistInduct = sigmaEffResistInduct;
        lumped.currentCoeff = currentCoeff;

        double G1_usual =
            (1.0 - sigma * dt / (2.0 * epsilon)) / (1.0 + sigma * dt / (2.0 * epsilon));
        double G2_usual =
            dt / epsilon / (1.0 + sigma * dt / (2.0 * epsilon));
        if (G1_usual < 0.0) {
            G1_usual = std::exp(-sigma * dt / epsilon);
            G2_usual = (1.0 - G1_usual) / sigma;
        }
        lumped.G1_usual = G1_usual;
        lumped.G2a_usual = G2_usual / transversalDeltaHa;
        lumped.G2b_usual = G2_usual / transversalDeltaHb;
    }
}

void LumpedSolver_t::advance(int timestep, double dt) {
    for (auto& lumped : nodes) {
        if (!lumped.Efield) {
            continue;
        }
        const auto& m = lumped.mat;

        if (m.inductor) {
            lumped.EfieldPrevPrev = lumped.EfieldPrev;
            lumped.EfieldPrev = *lumped.Efield;
            lumped.Jcur = lumped.currentCoeff * lumped.Jcur +
                          lumped.sigmaEffResistInduct *
                              (lumped.EfieldPrev + lumped.EfieldPrevPrev);
        } else {
            lumped.Jcur = 0.0;
        }

        const double dHa = *lumped.Ha_Plus - *lumped.Ha_Minu;
        const double dHb = *lumped.Hb_Plus - *lumped.Hb_Minu;
        const double curlH = lumped.G2a * dHa - lumped.G2b * dHb;

        if (m.resistor) {
            const double time = timestep * dt;
            if (time >= m.Rtime_on && time <= m.Rtime_off) {
                *lumped.Efield =
                    lumped.G1 * *lumped.Efield + curlH - lumped.GJ * lumped.Jcur;
            } else {
                *lumped.Efield = lumped.G1_usual * *lumped.Efield +
                                 lumped.G2a_usual * dHa - lumped.G2b_usual * dHb;
            }
        } else if (!m.diodo) {
            *lumped.Efield = lumped.G1 * *lumped.Efield + curlH - lumped.GJ * lumped.Jcur;
        }
    }
}

} // namespace lumped_slim
