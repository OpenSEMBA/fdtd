#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <stdexcept>

// Assuming these types and functions are defined elsewhere based on the Fortran context
// FDETYPES_m
using RKIND = double;
const RKIND RKIND_VAL = 1.0; // Placeholder for kind parameter usage

// Report_m
void StopOnError(int code1, int code2, const std::string& message) {
    std::cerr << "Error: " << message << std::endl;
    // In a real scenario, this might throw or exit
    throw std::runtime_error(message);
}

// Derived types from Fortran (simplified representation)
struct MediaInfo {
    RKIND SigmaM;
    RKIND Epr;
    RKIND Mur;
    RKIND Sigma;
    struct {
        bool already_YEEadvanced_byconformal;
        bool split_and_useless;
        bool ConformalPEC;
        bool multiport;
        bool AnisMultiport;
        bool pec;
        bool lumped;
        bool SGBC;
        bool Anisotropic;
        bool EDispersive;
        bool MDispersive;
        bool EdispersiveANIS;
        bool MdispersiveANIS;
    } Is;
    // Placeholder for multiport data if needed, though not fully expanded in this snippet
    struct {
        int numcapas;
        std::vector<RKIND> width;
        std::vector<RKIND> sigma;
        std::vector<RKIND> epr;
    } multiport;
};

struct SGGFDTDINFO_t {
    int NumMedia;
    RKIND dt;
    std::vector<MediaInfo> Med;
};

struct constants_t {
    std::vector<RKIND> g1;
    std::vector<RKIND> g2;
    std::vector<RKIND> gm1;
    std::vector<RKIND> gm2;
};

namespace CALC_CONSTANTS_m {

    void calc_g1g2gm1gm2(const SGGFDTDINFO_t& sgg, constants_t& g, RKIND& Eps0, RKIND& Mu0) {
        int r, i;
        RKIND Sigmam, Epsilon, Mu, Sigma, width, epr;
        // Character buffer not strictly needed for logic, but kept for potential error reporting if needed
        // std::string buff; 

        // Ensure vectors are sized correctly
        if (g.g1.size() != sgg.NumMedia + 1) {
            g.g1.resize(sgg.NumMedia + 1);
            g.g2.resize(sgg.NumMedia + 1);
            g.gm1.resize(sgg.NumMedia + 1);
            g.gm2.resize(sgg.NumMedia + 1);
        }

        for (r = 0; r <= sgg.NumMedia; ++r) {
            Sigmam = sgg.Med[r].SigmaM;
            Epsilon = Eps0 * sgg.Med[r].Epr;
            Mu = Mu0 * sgg.Med[r].Mur;
            Sigma = sgg.Med[r].Sigma;

            // In case of Multiport set updating Coefficients trivially, since the field
            // updating are handled locally in the Multiport module
            if ((sgg.Med[r].Is.already_YEEadvanced_byconformal) || (sgg.Med[r].Is.split_and_useless)) {
                g.g1[r] = 1.0;
                g.g2[r] = 0.0;
                g.gm1[r] = 1.0;
                g.gm2[r] = 0.0;
            } else {
                if (sgg.Med[r].Is.ConformalPEC) {
                    g.g1[r] = 1.0;
                    g.g2[r] = sgg.dt / Epsilon;
                    g.gm1[r] = 1.0;
                    g.gm2[r] = sgg.dt / Mu;
                } else if ((sgg.Med[r].Is.multiport) || (sgg.Med[r].Is.AnisMultiport)) {
                    g.g1[r] = 0.0; // null fields on the main procedure (good both for Ian and for me)
                    g.g2[r] = 0.0;
                    g.gm1[r] = 0.0;
                    g.gm2[r] = 0.0;
                } else if ((sgg.Med[r].Is.pec) || (r == 0)) {
                    // Trivially PEC updating Coefficients are set to 0.0
                    g.g1[r] = 0.0;
                    g.g2[r] = 0.0;
                    g.gm1[r] = 0.0;
                    g.gm2[r] = 0.0;
                } else if (sgg.Med[r].Is.lumped) {
                    // Trivially PEC NOT updating coefficients para el avance en E. They are set to 1.0. La rutina propia ya se encargara
                    g.g1[r] = 1.0;
                    g.g2[r] = 0.0;
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    if (g.gm1[r] < 0.0) { // exponential time stepping
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                } else if (sgg.Med[r].Is.SGBC) {
                    // 090519 He quitado todo este calculo que luego hara InitSGBCs para no duplicar codigo propenso a errores. Uso valores absurdos por lo que truene.
                    // ojo que los parametros stochastic tambien se obtendran en InitSGBCs, por eso lo he quitado esto de aqui
                    g.g1[r] = 0.0;
                    g.g2[r] = 0.0;
                    // g%g1(r)=1.0e31; g%g2(r)=-1.0e23 !quitado este sinsentido por si diese lugar a errores de redondeo 170519
                    
                    // quitado soporte multicapas en SGBC a 260419
                    // los intrinsecos suyos promediados si es una multicapa!!. La rutina de avance machaca los campos pero le paso estos coefs a initSGBCs para que los aproveche
                    
                    // hasta aqui lo comentado a 090519. Los calculos de g%gm1 y g%gm2 no los hace InitSGBCs y sus incertidumbres son nulas. asi que los hago aqui
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    if (g.gm1[r] < 0.0) { // exponential time stepping
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                } else if (sgg.Med[r].Is.Anisotropic) {
                    g.g1[r] = 1.0; // para que no haga nada en el bucle principal evitando los ifs
                    g.g2[r] = 0.0;
                    g.gm1[r] = 1.0; // para que no haga nada en el bucle principal evitando los ifs
                    g.gm2[r] = 0.0;
                } else if ((sgg.Med[r].Is.EDispersive) && (!sgg.Med[r].Is.MDispersive) && (!sgg.Med[r].Is.EdispersiveANIS)) {
                    // solo cierto para ISOTROPOS
                    g.g1[r] = 0.0; // will be overwritten by own values created by InitEDispersives
                    g.g2[r] = 0.0; // will be overwritten by own values created by InitEDispersives
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    if (g.gm1[r] < 0.0) { // exponential time stepping
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                } else if ((sgg.Med[r].Is.MDispersive) && (!sgg.Med[r].Is.EDispersive) && (!sgg.Med[r].Is.MdispersiveANIS)) {
                    // solo cierto para ISOTROPOS
                    g.g1[r] = (1.0 - Sigma * sgg.dt / (2.0 * Epsilon)) / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    g.g2[r] = sgg.dt / Epsilon / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    if (g.g1[r] < 0.0) { // exponential time stepping
                        g.g1[r] = std::exp(-Sigma * sgg.dt / Epsilon);
                        g.g2[r] = (1.0 - g.g1[r]) / Sigma;
                    }
                    g.gm1[r] = 0.0; // will be overwritten by own values created by InitMDispersives !when written this routine
                    g.gm2[r] = 0.0; // will be overwritten by own values created by InitMDispersives
                } else if ((sgg.Med[r].Is.MdispersiveANIS) || (sgg.Med[r].Is.EdispersiveANIS)) {
                    std::string buff = "ERROR: ANISOTROPIC DISPERSIVE CURRENTLY UNSUPPORTED IN THE ENGINE";
                    StopOnError(0, 0, buff); // lo deberia reportar y parar antes SEMBA_FDTD.F90 !quitar algun dia para que no ralentice 170719
                } else {
                    g.g1[r] = (1.0 - Sigma * sgg.dt / (2.0 * Epsilon)) / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    g.g2[r] = sgg.dt / Epsilon / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    if (g.g1[r] < 0.0) { // exponential time stepping
                        g.g1[r] = std::exp(-Sigma * sgg.dt / Epsilon);
                        g.g2[r] = (1.0 - g.g1[r]) / Sigma;
                    }
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    if (g.gm1[r] < 0.0) { // exponential time stepping
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                }
            }
        }
    }

} // namespace CALC_CONSTANTS_m