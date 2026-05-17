#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <stdexcept>

// Forward declarations for types defined in other modules
// These would typically be in headers like FDETYPES_m.h and Report_m.h

// Placeholder for FDETYPES_m types
struct SGGFDTDINFO_t;
struct constants_t;

// Placeholder for Report_m functions
void StopOnError(int code1, int code2, const std::string& message);

// Placeholder for constants from FDETYPES_m
const double RKIND = 1.0; // Assuming RKIND maps to double
const int BUFSIZE = 256;  // Assuming BUFSIZE is a constant buffer size

namespace CALC_CONSTANTS_m {

    void calc_g1g2gm1gm2(const SGGFDTDINFO_t& sgg, constants_t& g, double& Eps0, double& Mu0) {
        int r, i;
        double Sigmam, Epsilon, Mu, Sigma, width, epr;
        std::string buff;

        // Fortran do loops are inclusive. Assuming sgg.NumMedia is the upper bound.
        // Note: Fortran arrays are often 1-based or 0-based depending on declaration.
        // The code uses r=0 to sgg%NumMedia. We assume the vector/array access handles this or
        // that the underlying data structure supports this indexing.
        // If sgg.Med is a vector, we must ensure it has size at least NumMedia + 1.
        
        for (r = 0; r <= sgg.NumMedia; ++r) {
            // Accessing members. Assuming sgg.Med is a vector or array.
            // In C++, if Med is std::vector, we use sgg.Med[r].
            // If it's a raw array, we use sgg.Med[r].
            // We assume the struct SGGFDTDINFO_t has a member 'Med' that is indexable.
            
            // Note: The original code uses sgg%Med(r). In C++, this translates to sgg.Med[r].
            // We assume sgg.Med is a std::vector<Media_t> or similar.
            
            Sigmam = sgg.Med[r].SigmaM;
            Epsilon = Eps0 * sgg.Med[r].Epr;
            Mu = Mu0 * sgg.Med[r].Mur;
            Sigma = sgg.Med[r].Sigma;

            // Check conditions. Note: Fortran logicals map to bool.
            // The original code uses sgg%Med(R)%Is%...
            // We assume sgg.Med[r].Is is a struct containing boolean flags.
            
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
                    g.g1[r] = 0.0;
                    g.g2[r] = 0.0;
                    g.gm1[r] = 0.0;
                    g.gm2[r] = 0.0;
                } else if ((sgg.Med[r].Is.pec) || (r == 0)) {
                    g.g1[r] = 0.0;
                    g.g2[r] = 0.0;
                    g.gm1[r] = 0.0;
                    g.gm2[r] = 0.0;
                } else if (sgg.Med[r].Is.lumped) {
                    g.g1[r] = 1.0;
                    g.g2[r] = 0.0;
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    
                    if (g.gm1[r] < 0.0) {
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                } else if (sgg.Med[r].Is.SGBC) {
                    // Commented out code removed.
                    g.g1[r] = 0.0;
                    g.g2[r] = 0.0;
                    
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    
                    if (g.gm1[r] < 0.0) {
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                } else if (sgg.Med[r].Is.Anisotropic) {
                    g.g1[r] = 1.0;
                    g.g2[r] = 0.0;
                    g.gm1[r] = 1.0;
                    g.gm2[r] = 0.0;
                } else if ((sgg.Med[r].Is.EDispersive) && (!sgg.Med[r].Is.MDispersive) && (!sgg.Med[r].Is.EdispersiveANIS)) {
                    g.g1[r] = 0.0;
                    g.g2[r] = 0.0;
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    
                    if (g.gm1[r] < 0.0) {
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                } else if ((sgg.Med[r].Is.MDispersive) && (!sgg.Med[r].Is.EDispersive) && (!sgg.Med[r].Is.MdispersiveANIS)) {
                    g.g1[r] = (1.0 - Sigma * sgg.dt / (2.0 * Epsilon)) / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    g.g2[r] = sgg.dt / Epsilon / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    
                    if (g.g1[r] < 0.0) {
                        g.g1[r] = std::exp(-Sigma * sgg.dt / Epsilon);
                        g.g2[r] = (1.0 - g.g1[r]) / Sigma;
                    }
                    
                    g.gm1[r] = 0.0;
                    g.gm2[r] = 0.0;
                } else if ((sgg.Med[r].Is.MdispersiveANIS) || (sgg.Med[r].Is.EdispersiveANIS)) {
                    buff = "ERROR: ANISOTROPIC DISPERSIVE CURRENTLY UNSUPPORTED IN THE ENGINE";
                    StopOnError(0, 0, buff);
                } else {
                    g.g1[r] = (1.0 - Sigma * sgg.dt / (2.0 * Epsilon)) / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    g.g2[r] = sgg.dt / Epsilon / (1.0 + Sigma * sgg.dt / (2.0 * Epsilon));
                    
                    if (g.g1[r] < 0.0) {
                        g.g1[r] = std::exp(-Sigma * sgg.dt / Epsilon);
                        g.g2[r] = (1.0 - g.g1[r]) / Sigma;
                    }
                    
                    g.gm1[r] = (1.0 - Sigmam * sgg.dt / (2.0 * Mu)) / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    g.gm2[r] = sgg.dt / Mu / (1.0 + Sigmam * sgg.dt / (2.0 * Mu));
                    
                    if (g.gm1[r] < 0.0) {
                        g.gm1[r] = std::exp(-Sigmam * sgg.dt / Mu);
                        g.gm2[r] = (1.0 - g.gm1[r]) / Sigmam;
                    }
                }
            }
        }
    }

} // namespace CALC_CONSTANTS_m