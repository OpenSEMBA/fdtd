#ifndef OBSERVATION_TYPES_H
#define OBSERVATION_TYPES_H

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdint>

namespace Observa_m {

using RKIND = double;
using RKIND_tiempo = double;

struct Serialized_t {
    std::vector<double> Valor;
    std::vector<double> Valor_x;
    std::vector<double> Valor_y;
    std::vector<double> Valor_z;
    std::vector<double> ValorE;
    std::vector<double> Valor_Ex;
    std::vector<double> Valor_Ey;
    std::vector<double> Valor_Ez;
    std::vector<double> ValorH;
    std::vector<double> Valor_Hx;
    std::vector<double> Valor_Hy;
    std::vector<double> Valor_Hz;

    std::vector<std::complex<double>> ValorComplex_x;
    std::vector<std::complex<double>> ValorComplex_y;
    std::vector<std::complex<double>> ValorComplex_z;
    std::vector<std::complex<double>> ValorComplex_Ex;
    std::vector<std::complex<double>> ValorComplex_Ey;
    std::vector<std::complex<double>> ValorComplex_Ez;
    std::vector<std::complex<double>> ValorComplex_Hx;
    std::vector<std::complex<double>> ValorComplex_Hy;
    std::vector<std::complex<double>> ValorComplex_Hz;

    std::vector<int> eI;
    std::vector<int> eJ;
    std::vector<int> eK;
    std::vector<int> currentType;
    std::vector<int> sggmtag;

    void allocate_for_time_domain(int numberOfSerialized) {
        const size_t n = static_cast<size_t>(numberOfSerialized);
        Valor.assign(n, 0.0);
        Valor_x.assign(n, 0.0);
        Valor_y.assign(n, 0.0);
        Valor_z.assign(n, 0.0);
        ValorE.assign(n, 0.0);
        Valor_Ex.assign(n, 0.0);
        Valor_Ey.assign(n, 0.0);
        Valor_Ez.assign(n, 0.0);
        ValorH.assign(n, 0.0);
        Valor_Hx.assign(n, 0.0);
        Valor_Hy.assign(n, 0.0);
        Valor_Hz.assign(n, 0.0);
    }

    void deallocate_for_time_domain() {
        Valor.clear();
        Valor_x.clear();
        Valor_y.clear();
        Valor_z.clear();
        ValorE.clear();
        Valor_Ex.clear();
        Valor_Ey.clear();
        Valor_Ez.clear();
        ValorH.clear();
        Valor_Hx.clear();
        Valor_Hy.clear();
        Valor_Hz.clear();
    }

    void allocate_for_frequency_domain(int numberOfSerialized) {
        allocate_for_time_domain(numberOfSerialized);
        const size_t n = static_cast<size_t>(numberOfSerialized);
        ValorComplex_x.assign(n, {0.0, 0.0});
        ValorComplex_y.assign(n, {0.0, 0.0});
        ValorComplex_z.assign(n, {0.0, 0.0});
        ValorComplex_Ex.assign(n, {0.0, 0.0});
        ValorComplex_Ey.assign(n, {0.0, 0.0});
        ValorComplex_Ez.assign(n, {0.0, 0.0});
        ValorComplex_Hx.assign(n, {0.0, 0.0});
        ValorComplex_Hy.assign(n, {0.0, 0.0});
        ValorComplex_Hz.assign(n, {0.0, 0.0});
    }

    void deallocate_for_frequency_domain() {
        deallocate_for_time_domain();
        ValorComplex_x.clear();
        ValorComplex_y.clear();
        ValorComplex_z.clear();
        ValorComplex_Ex.clear();
        ValorComplex_Ey.clear();
        ValorComplex_Ez.clear();
        ValorComplex_Hx.clear();
        ValorComplex_Hy.clear();
        ValorComplex_Hz.clear();
    }

    void allocate_current_value(int numberOfSerialized) {
        const size_t n = static_cast<size_t>(numberOfSerialized);
        eI.assign(n, 0);
        eJ.assign(n, 0);
        eK.assign(n, 0);
        currentType.assign(n, 0);
        sggmtag.assign(n, 0);
    }

    void deallocate_current_value() {
        eI.clear();
        eJ.clear();
        eK.clear();
        currentType.clear();
        sggmtag.clear();
    }
};

struct output_t {
    bool SaveAll = false;
    int Trancos = 0;
};

struct Obses_t {
    RKIND_tiempo TimeStep = 0.0;
    RKIND_tiempo InitialTime = 0.0;
    RKIND_tiempo FinalTime = 0.0;
    int32_t nP = 0;
    bool Volumic = false;
    RKIND InitialFreq = 0.0;
    RKIND FinalFreq = 0.0;
    RKIND FreqStep = 0.0;
};

inline bool approx_equal(RKIND a, RKIND b, RKIND tol) {
    return std::abs(a - b) <= tol;
}

inline void preprocess_observation(Obses_t& obs, output_t& out,
                                   const std::vector<RKIND_tiempo>& /*time*/,
                                   int finaltimestep, RKIND_tiempo dt, bool saveall) {
    if (obs.InitialTime < obs.TimeStep) {
        obs.InitialTime = 0.0;
    }
    if (obs.TimeStep > dt) {
        obs.TimeStep = dt;
    }
    if (obs.FreqStep <= 0.0 || obs.FreqStep > (obs.FinalFreq - obs.InitialFreq)) {
        obs.FreqStep = obs.FinalFreq - obs.InitialFreq;
        if (obs.FreqStep <= 0.0) obs.FreqStep = 1.0;
    }
    if (!obs.Volumic && saveall) {
        out.SaveAll = true;
    }
    if (obs.FinalTime < obs.InitialTime) {
        obs.FinalTime = obs.InitialTime;
    }
    (void)finaltimestep;
}

} // namespace Observa_m

#endif
