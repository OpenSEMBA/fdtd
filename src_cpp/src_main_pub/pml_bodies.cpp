#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <cstring>

struct sim_control_t {
    int layoutnumber;
    int num_procs;
    bool resume;
};

struct SGGFDTDINFO_t {
    struct { int XI, XE, YI, YE, ZI, ZE; } alloc[10];
    struct { int XI, XE, YI, YE, ZI, ZE; } SINPMLSweep[10];
    int nummedia;
    struct {
        struct { bool PMLbody; struct { int orient; } PMLbody_arr[1]; } Is;
    } Med[100];
    double dt;
};

struct media_matrices_t {
    int sggMiEx[100][100][100];
    int sggMiEy[100][100][100];
    int sggMiEz[100][100][100];
    int sggMiHx[100][100][100];
    int sggMiHy[100][100][100];
    int sggMiHz[100][100][100];
};

typedef double RKIND;
const int iEx = 0, iEy = 1, iEz = 2;
const int iHx = 3, iHy = 4, iHz = 5;
const int BUFSIZE = 256;

void WarnErrReport(const std::string&, bool);
void print11(int, const std::string&);
const std::string SEPARADOR = "========================================";

namespace PMLbodies_m {

    struct BerPML__t {
        RKIND* field;
        RKIND* Plus;
        RKIND* Minu;
        RKIND gx2;
        RKIND P_be;
        RKIND P_ce;
        RKIND Psi;
        RKIND transversalDelta;
        RKIND del;
        RKIND posi;
        int minTotal;
        int maxTotal;
    };

    struct berpml_t {
        int NumNodes;
        int orient;
        std::vector<BerPML__t> nodes;
    };

    berpml_t berpmlE;
    berpml_t berpmlH;
    const int PMLorden = 2;
    const RKIND CoeffReflPML = 1e-4;
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    void InitPMLbodies(
        const SGGFDTDINFO_t&, const media_matrices_t&,
        const RKIND*, const RKIND*, const RKIND*,
        const RKIND*, const RKIND*, const RKIND*,
        const int*, const int*, const int*,
        const int*, const int*, const int*,
        const RKIND*, const RKIND*,
        bool&, const sim_control_t&, RKIND, RKIND) {}
    void calc_pmlbodypar(const SGGFDTDINFO_t&, RKIND, RKIND) {}
    void AdvancePMLbodyE() {}
    void AdvancePMLbodyH() {}
    void StorefieldsPMLbodies() {}
    void DestroyPMLbodies(SGGFDTDINFO_t&) {}

}
