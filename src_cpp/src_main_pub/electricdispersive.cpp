#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>

struct SGGFDTDINFO_t;
struct media_matrices_t;

void print11(int, const std::string&) {}
const std::string SEPARADOR = "========================================";

namespace EDispersives_m {

    using RKIND = double;
    using CKIND = std::complex<double>;

    extern const int iEx;
    extern const int iEy;
    extern const int iEz;

    struct field_t {
        int32_t i;
        int32_t j;
        int32_t k;
        int32_t WhatField;
        RKIND* FieldPresent;
        RKIND FieldPrevious;
        std::vector<CKIND> Current;
    };

    struct EDispersive_t {
        int32_t indexmed;
        int32_t numnodesEx;
        int32_t numnodesEy;
        int32_t numnodesEz;
        int32_t numpolres11;
        std::vector<CKIND> Beta;
        std::vector<CKIND> Kappa;
        std::vector<CKIND> G3;
        std::vector<field_t> NodesEx;
        std::vector<field_t> NodesEy;
        std::vector<field_t> NodesEz;
    };

    struct EDispersive2_t {
        int32_t NumEDispersives;
        std::vector<EDispersive_t> Medium;
    };

    EDispersive2_t Dutton;

    void InitEDispersives(const SGGFDTDINFO_t&, const media_matrices_t&, std::vector<RKIND>&, std::vector<RKIND>&, std::vector<std::vector<std::vector<RKIND>>>&, std::vector<std::vector<std::vector<RKIND>>>&, std::vector<std::vector<std::vector<RKIND>>>&, bool&, bool) {}
    void AdvanceEDispersiveE(const SGGFDTDINFO_t&) {}
    void StoreFieldsEDispersives() {}
    void DestroyEDispersives(SGGFDTDINFO_t&) {}

}
