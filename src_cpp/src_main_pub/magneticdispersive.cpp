#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>

struct SGGFDTDINFO_t;
struct media_matrices_t;

namespace Mdispersives_m {

    struct field_t {
        int32_t i, j, k;
        int32_t WhatField;
        double* FieldPresent;
        double FieldPrevious;
        std::vector<std::complex<double>> Current;
    };

    struct Mdispersive_t {
        int32_t indexmed;
        int32_t numnodesHx, numnodesHy, numnodesHz;
        int32_t numpolres11;
        std::vector<std::complex<double>> Beta;
        std::vector<std::complex<double>> Kappa;
        std::vector<std::complex<double>> GM3;
        std::vector<field_t> NodesHx;
        std::vector<field_t> NodesHy;
        std::vector<field_t> NodesHz;
    };

    struct Mdispersive2_t {
        int32_t NumMdispersives;
        std::vector<Mdispersive_t> Medium;
    };

    Mdispersive2_t MDutton;

    void InitMdispersives(SGGFDTDINFO_t&, media_matrices_t&, std::vector<double>&, std::vector<double>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, std::vector<std::vector<std::vector<double>>>&, bool&, bool) {}
    void AdvanceMdispersiveH(SGGFDTDINFO_t&) {}
    void StoreFieldsMdispersives() {}
    void DestroyMdispersives(SGGFDTDINFO_t&) {}

}
